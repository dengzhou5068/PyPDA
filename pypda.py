import argparse
import os
import requests
import configparser
import json
import warnings
import shutil
from pathlib import Path
import xml.etree.ElementTree as ET
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio import ExPASy
from Bio import SwissProt
from Bio.PDB.PDBList import PDBList
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from glob import glob
from concurrent.futures import ThreadPoolExecutor
from functools import lru_cache

# 过滤 PDB 构建警告
warnings.filterwarnings("ignore", category=PDBConstructionWarning)

# 读取序列配置文件
seq_config = configparser.ConfigParser()
seq_config.read('seq_config.ini')

# 读取PDB配置文件
pdb_config = configparser.ConfigParser()
pdb_config.read('pdb_config.ini')

# 从配置文件中获取API URL
UNIPROT_API = seq_config.get('UNIPROT_API', 'UNIPROT_API')
UNIPROT_XML_API = seq_config.get('UNIPROT_API', 'UNIPROT_XML_API')

# 从配置文件中获取命名空间映射
ns = dict(seq_config.items('XML_NAMESPACES'))
ERROR_LOG_PATH = os.path.join(os.getcwd(), 'error.txt')

# ------------------------------ common_utils 功能 ------------------------------

def fetch_protein_sequences(genes, output_dir):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    for gene in genes:
        params = {
            "query": f"gene_exact:{gene} AND organism_id:9606",
            "format": "fasta",
            "fields": "sequence"
        }
        try:
            response = requests.get(UNIPROT_API, params=params)
            response.raise_for_status()
            if not response.text:
                print(f"未找到基因 {gene} 对应的蛋白质序列")
                continue
            output_file = output_dir / f"{gene}.fasta"
            with open(output_file, "w", encoding="utf-8") as f:
                first_seq = response.text.split('>', 2)[1] if '>' in response.text else response.text
                f.write(f">{first_seq}")
            print(f"已将 {gene} 的第一个全长蛋白质序列保存到 {output_file}")
            uniprot_id = first_seq.split('\n', 1)[0].split('|')[1]
            yield gene, uniprot_id
        except requests.exceptions.RequestException as e:
            print(f"获取 {gene} 失败: {str(e)}")

def fetch_domain_information(gene_uniprot_pairs, output_dir):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    domain_info_file = output_dir / "domain_info.md"
    with open(domain_info_file, "w", encoding="utf-8") as domain_f:
        domain_f.write("# 结构域信息\n\n")
    for gene, uniprot_id in gene_uniprot_pairs:
        try:
            fasta_file = output_dir / f"{gene}.fasta"
            if fasta_file.exists():
                record = SeqIO.read(fasta_file, "fasta")
                total_amino_acids = len(record.seq)
            else:
                total_amino_acids = "未知"
            
            xml_response = requests.get(f"{UNIPROT_XML_API}{uniprot_id}.xml")
            xml_response.raise_for_status()
            root = ET.fromstring(xml_response.text)
            with open(domain_info_file, "a", encoding="utf-8") as domain_f:
                domain_f.write(f"## {gene}\n")
                domain_f.write(f"### 总氨基酸数量: {total_amino_acids}\n\n")
                for feature in root.findall('.//uniprot:feature[@type="domain"]', ns):
                    begin = int(feature.find('uniprot:location/uniprot:begin', ns).attrib['position'])
                    end = int(feature.find('uniprot:location/uniprot:end', ns).attrib['position'])
                    domain_name = feature.attrib.get('description', '未知结构域')
                    domain_f.write(f"- 结构域名称: {domain_name}, 序列范围: {begin}-{end}\n")
                    print(f"{gene} 的 {domain_name} 结构域的序列编号范围: {begin}-{end}")
        except requests.exceptions.RequestException as e:
            print(f"获取 {gene} 的结构域信息失败: {str(e)}")

def read_fasta(file_path, return_header=False):
    """读取FASTA文件，可选择返回头部信息
    
    Args:
        file_path: 文件路径
        return_header: 是否返回头部信息
    
    Returns:
        序列字符串，或(头部, 序列)元组
    """
    record = SeqIO.read(file_path, "fasta")
    if return_header:
        return record.description, str(record.seq)
    return str(record.seq)

def parse_arguments(description, args_config):
    parser = argparse.ArgumentParser(description=description)
    for arg in args_config:
        if 'nargs' in arg:
            parser.add_argument(arg['name'], nargs=arg['nargs'], type=arg.get('type'), required=arg.get('required', False), help=arg['help'])
        else:
            parser.add_argument(arg['name'], type=arg.get('type'), required=arg.get('required', False), help=arg['help'])
    return parser.parse_args()

def read_fasta_file(fasta_file):
    sequence = ""
    with open(fasta_file, 'r') as file:
        lines = file.readlines()
        header = lines[0].strip()
        for line in lines[1:]:
            sequence += line.strip()
    return header, sequence

def generate_output_filename(base_name, ext, *args):
    """统一生成输出文件名
    
    Args:
        base_name: 基础文件名
        ext: 文件扩展名（含.）
        *args: 可变参数，用于生成附加信息
    """
    if not args:
        return f"{base_name}{ext}"
    info = "_".join(map(str, args)) if len(args) > 1 else str(args[0])
    return f"{base_name}_{info}{ext}"

def save_sequence(header, sequence, output_file):
    with open(output_file, 'w', encoding="utf-8") as file:
        file.write(header + '\n')
        for i in range(0, len(sequence), 60):
            file.write(sequence[i:i+60] + '\n')

def validate_input(mutation_positions, new_amino_acids):
    if len(mutation_positions) != len(new_amino_acids):
        print("The number of mutation positions and new amino acids must be the same.")
        exit(1)
    return True

def perform_mutations(record, mutation_positions, new_amino_acids):
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    sequence = str(record.seq)
    for pos, aa in zip(mutation_positions, new_amino_acids):
        sequence = sequence[:pos - 1] + aa + sequence[pos:]
    return SeqRecord(Seq(sequence), id=record.id, description=record.description)

def extract_subsequence(fasta_file, start, end):
    header, sequence = read_fasta_file(fasta_file)
    if start < 1 or end > len(sequence):
        print('提取位置超出序列范围')
        return None
    subsequence = sequence[start - 1:end]
    return subsequence

# ------------------------------ seq 功能 ------------------------------

def generate_output_filename(fasta_file, *args):
    """统一生成输出文件名"""
    base_name = os.path.splitext(fasta_file)[0]
    if len(args) == 2:  # extract命令
        return f"{base_name}_{args[0]}-{args[1]}.fasta"
    else:  # mutate命令
        mutation_info = "".join([f"{pos}{aa}" for pos, aa in zip(args[0], args[1])])
        return f"{base_name}_{mutation_info}.fasta"

def read_align_fasta(file_path):
    """读取比对用的FASTA文件"""
    return str(SeqIO.read(file_path, "fasta").seq)

def pairwise_alignment(seq1, seq2):
    """执行双序列比对并计算同源性"""
    aligner = PairwiseAligner()
    best_alignment = next(aligner.align(seq1, seq2))
    score = best_alignment.score
    homology = (score / max(len(seq1), len(seq2))) * 100
    return score, homology

def compare_sequences(file_paths):
    """比较多个FASTA文件中的序列"""
    sequences = [read_align_fasta(fp) for fp in file_paths]
    for i, seq1 in enumerate(sequences):
        for j, seq2 in enumerate(sequences[i+1:], i+1):
            score, homology = pairwise_alignment(seq1, seq2)
            print(f"比对 {i+1} 和 {j+1}:\n比对得分: {score}\n同源性百分比: {homology:.2f}%\n")

def seq_main(args):
    if args.command == "fetch":
        with open(args.input_file) as f:
            genes = [line.strip() for line in f if line.strip()]
        if genes:
            gene_uniprot_pairs = list(fetch_protein_sequences(genes, args.output_dir))
            fetch_domain_information(gene_uniprot_pairs, args.output_dir)
        else:
            print("输入文件中未找到有效的基因名称")

    elif args.command == "extract":
        if not os.path.isfile(args.fasta_file):
            print(f"错误：文件 {args.fasta_file} 不存在")
            return
        header, _ = read_fasta_file(args.fasta_file)
        extracted_seq = extract_subsequence(args.fasta_file, args.start, args.end)
        output_file = generate_output_filename(args.fasta_file, args.start, args.end)
        save_sequence(header, extracted_seq, output_file)
        print(f"已保存截取序列至 {output_file}")

    elif args.command == "mutate":
        validate_input(args.pos, args.aa)
        record = next(SeqIO.parse(args.fasta_file, "fasta"))
        mutated_record = perform_mutations(record, args.pos, args.aa)
        output_file = generate_output_filename(args.fasta_file, args.pos, args.aa)
        SeqIO.write(mutated_record, output_file, "fasta")
        print(f"已保存突变序列至 {output_file}")

    elif args.command == "align":
        compare_sequences(args.fasta_files)

# ------------------------------ pdbquery 功能 ------------------------------

def log_error(error_msg, error_file_path):
    """
    公共错误记录函数，用于统一写入错误信息到指定文件。
    
    :param error_msg: 错误信息字符串
    :param error_file_path: 错误文件路径
    """
    with open(error_file_path, 'a', encoding='utf-8') as f:
        f.write(f"{error_msg}\n")
    print(error_msg)

def parse_exclude_residues():
    """
    解析排除残基的配置文件，获取需要排除的残基列表。

    :return: 需要排除的残基列表
    """
    config = configparser.ConfigParser()
    config.read('exclude_residues.ini')
    exclude_residues = []
    for section in config.sections():
        for key in config[section]:
            exclude_residues.extend([res.strip() for res in config[section][key].split(',')])
    return exclude_residues

@lru_cache(maxsize=None)
def get_pdb_ids_from_uniprot(uniprot_id):
    """
    通过UniProt ID获取对应的PDB ID列表（带缓存优化）。

    :param uniprot_id: UniProt ID
    :return: PDB ID列表
    """
    try:
        handle = ExPASy.get_sprot_raw(uniprot_id)
        record = SwissProt.read(handle)
        pdb_ids = []
        for cross_ref in record.cross_references:
            if cross_ref[0] == 'PDB':
                pdb_ids.append(cross_ref[1])
        return pdb_ids
    except Exception as e:
        log_error(f"获取PDB ID时出错: {e}", ERROR_LOG_PATH)
        return []

# 添加在文件顶部工具函数区域

def parallel_executor(func, items, max_workers=None):
    """并行执行函数
    
    Args:
        func: 要并行执行的函数
        items: 迭代参数列表
        max_workers: 最大工作线程数
    """
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        executor.map(func, items)

def download_pdb_files(pdb_ids, output_dir):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    pdbl = PDBList()
    failed_files = []
    error_file = output_dir / 'error.txt'

    def download_single_file(pdb_id):
        try:
            pdbl.retrieve_pdb_file(pdb_id, pdir=output_dir, file_format='mmCif')
            print(f"成功下载CIF文件: {pdb_id}.cif")
        except Exception as e:
            error_msg = f"下载CIF文件 {pdb_id} 时出错: {e}"
            print(error_msg)
            failed_files.append(pdb_id)
            with open(error_file, 'a', encoding='utf-8') as f:
                f.write(error_msg + '\n')

    # 修改download_pdb_files中的调用
    parallel_executor(download_single_file, pdb_ids)

    if failed_files:
        with open(os.path.join(output_dir, 'download_failed.txt'), 'a', encoding='utf-8') as f:
            for failed_file in failed_files:
                f.write(f"{failed_file}.cif\n")

def extract_ligands_from_pdb(folder_path, exclude_residues):
    """
    从.cif文件中提取配体信息。

    :param folder_path: 文件夹路径
    :param exclude_residues: 需要排除的残基列表
    :return: 配体到PDB ID的映射字典
    """
    parser = MMCIFParser()
    ligand_pdb_dict = {}
    error_file = os.path.join(folder_path, 'error.txt')

    for filename in os.listdir(folder_path):
        if filename.endswith('.cif'):
            pdb_id = filename.split('.')[0]
            try:
                structure = parser.get_structure(pdb_id, os.path.join(folder_path, filename))
                ligands = []
                for model in structure:
                    for chain in model:
                        for residue in chain:
                            if residue.id[0].startswith('H_'):
                                residue_name = residue.resname.strip()
                                if residue_name not in exclude_residues and residue_name not in ligands:
                                    ligands.append(residue_name)
            except Exception as e:
                error_msg = f"解析 {filename} 时出错: {e}"
                print(error_msg)
                with open(error_file, 'a', encoding='utf-8') as f:
                    f.write(error_msg + '\n')
                continue

            for ligand in ligands:
                if ligand in ligand_pdb_dict:
                    ligand_pdb_dict[ligand].append(pdb_id)
                else:
                    ligand_pdb_dict[ligand] = [pdb_id]
    return ligand_pdb_dict

def write_ligand_info_to_md(ligand_pdb_dict, output_dir):
    """
    将配体信息写入pdb_ligand.md文件。

    :param ligand_pdb_dict: 配体到PDB ID的映射字典
    :param output_dir: 输出目录
    """
    output_file = 'pdb_ligand.md'
    with open(os.path.join(output_dir, output_file), 'w') as outfile:
        outfile.write('| Ligands | PDB ID |\n')
        outfile.write('| --- | --- |\n')
        for ligand, pdb_ids in ligand_pdb_dict.items():
            pdb_id_str = ', '.join(pdb_ids)
            outfile.write(f'| {ligand} | {pdb_id_str} |\n')
    print(f'结果已写入 {os.path.join(output_dir, output_file)}')

def move_no_ligand_files(folder_path, ligand_pdb_dict):
    """
    将没有配体的CIF文件移动到no_ligand文件夹。
    """
    no_ligand_dir = os.path.join(folder_path, 'no_ligand')
    if not os.path.exists(no_ligand_dir):
        os.makedirs(no_ligand_dir)

    pdb_ids_in_md = set()
    for pdb_ids in ligand_pdb_dict.values():
        pdb_ids_in_md.update(pdb_ids)

    error_file = os.path.join(folder_path, 'error.txt')
    for filename in os.listdir(folder_path):
        if filename.endswith('.cif'):
            pdb_id = filename.split('.')[0]
            if pdb_id not in pdb_ids_in_md:
                src_path = os.path.join(folder_path, filename)
                dst_path = os.path.join(no_ligand_dir, filename)
                try:
                    shutil.move(src_path, dst_path)
                    print(f"已将 {filename} 移动到 {no_ligand_dir}")
                except Exception as e:
                    log_error(f"移动 {filename} 时出错: {e}", error_file)

def move_with_ligand_files(folder_path, ligand_pdb_dict):
    """
    将含有配体的CIF文件移动到with_ligands文件夹。
    """
    with_ligands_dir = os.path.join(folder_path, 'with_ligands')
    if not os.path.exists(with_ligands_dir):
        os.makedirs(with_ligands_dir)

    pdb_ids_with_ligands = set()
    for pdb_ids in ligand_pdb_dict.values():
        pdb_ids_with_ligands.update(pdb_ids)

    error_file = os.path.join(folder_path, 'error.txt')
    for filename in os.listdir(folder_path):
        if filename.endswith('.cif'):
            pdb_id = filename.split('.')[0]
            if pdb_id in pdb_ids_with_ligands:
                src_path = os.path.join(folder_path, filename)
                dst_path = os.path.join(with_ligands_dir, filename)
                try:
                    shutil.move(src_path, dst_path)
                    print(f"已将含配体文件 {filename} 移动到 {with_ligands_dir}")
                except Exception as e:
                    log_error(f"移动含配体文件 {filename} 时出错: {e}", error_file)

def download_ligand_json(unique_ligands, output_dir):
    """
    通过API并行查询配体的JSON文件并下载保存到json子文件夹。
    """
    json_dir = os.path.join(output_dir, 'json')
    os.makedirs(json_dir, exist_ok=True)

    def download_single_ligand(ligand):
        retry = Retry(total=3, backoff_factor=1, status_forcelist=[500, 502, 503, 504])
        adapter = HTTPAdapter(max_retries=retry)
        s = requests.Session()
        s.mount('http://', adapter)
        s.mount('https://', adapter)
        
        try:
            url = f"https://data.rcsb.org/rest/v1/core/chemcomp/{ligand}"
            response = s.get(url, timeout=10)
            if response.status_code == 200:
                data = response.json()
                json_path = os.path.join(json_dir, f'{ligand}.json')
                with open(json_path, 'w', encoding='utf-8') as json_file:
                    json.dump(data, json_file, ensure_ascii=False, indent=4)
                print(f"JSON file for {ligand} has been saved to {json_path}")
            else:
                print(f"Request for {ligand} failed, status code: {response.status_code}")
        except Exception as e:
            print(f"下载配体 {ligand} 时出错: {e}")
        finally:
            s.close()

    with ThreadPoolExecutor() as executor:
        executor.map(download_single_ligand, unique_ligands)

def write_chemical_info_to_md(unique_ligands, output_dir):
    """
    从json子文件夹中读取JSON文件并写入chemical_components_info.md。
    """
    json_dir = os.path.join(output_dir, 'json')
    with open(os.path.join(output_dir, 'chemical_components_info.md'), 'w', encoding='utf-8') as md_file:
        md_file.write("| Chemical Component ID | name | formula | formula_weight | Canonical Smiles |\n")
        md_file.write("| --- | --- | --- | --- | --- |\n")
        for ligand in unique_ligands:
            try:
                json_path = os.path.join(json_dir, f'{ligand}.json')
                with open(json_path, 'r', encoding='utf-8') as json_file:
                    data = json.load(json_file)
                chem_comp = data.get('chem_comp', {})
                name = chem_comp.get('name', '')
                formula = chem_comp.get('formula', '')
                formula_weight = chem_comp.get('formula_weight', '')
                canonical_smiles = ''
                descriptors = data.get('pdbx_chem_comp_descriptor', [])
                for descriptor in descriptors:
                    if descriptor.get('type') == 'SMILES_CANONICAL':
                        canonical_smiles = descriptor.get('descriptor', '')
                        break
                md_file.write(f"| {ligand} | {name} | {formula} | {formula_weight} | {canonical_smiles} |\n")
                print(f"Information for {ligand} has been added to the MD file.")
            except FileNotFoundError:
                print(f"JSON file for {ligand} not found in {json_dir}")
    print(f"MD file {os.path.join(output_dir, 'chemical_components_info.md')} has been created successfully.")

def extract_ligands_coordinates(output_dir, exclude_residues):
    """
    直接通过文本解析.cif文件提取配体坐标并保存成cif文件（并行版）。
    """
    with_ligands_dir = os.path.join(output_dir, 'with_ligands')
    ligands_dir = os.path.join(with_ligands_dir, 'ligands')
    os.makedirs(ligands_dir, exist_ok=True)
    
    cif_files = glob(os.path.join(with_ligands_dir, '*.cif'))
    error_file = os.path.join(output_dir, 'error.txt')

    def parse_cif_atom_sites(cif_content):
        atom_site_fields = []
        atom_site_data = []
        in_loop = False
        
        lines = cif_content.split('\n')
        for line in lines:
            line = line.strip()
            if line.startswith('_atom_site.'):
                in_loop = True
                atom_site_fields.append(line.split('.')[1])
            elif in_loop:
                if line.startswith('#'):
                    in_loop = False
                elif line:
                    atom_site_data.append(line.split())
        
        field_index = {field: idx for idx, field in enumerate(atom_site_fields)}
        required_fields = ['group_PDB', 'auth_comp_id', 'auth_asym_id', 'auth_seq_id', 
                          'Cartn_x', 'Cartn_y', 'Cartn_z', 'id', 'type_symbol']
        
        if not all(f in field_index for f in required_fields):
            return None, "缺少必要的_atom_site字段"
        
        ligands = {}
        for data in atom_site_data:
            if data[field_index['group_PDB']] != 'HETATM':
                continue
            
            residue_name = data[field_index['auth_comp_id']].strip()
            if residue_name in exclude_residues:
                continue
            
            try:
                x = float(data[field_index['Cartn_x']])
                y = float(data[field_index['Cartn_y']])
                z = float(data[field_index['Cartn_z']])
            except (ValueError, IndexError):
                continue
            
            atom_fields = [
                data[field_index['group_PDB']],
                residue_name,
                data[field_index['auth_asym_id']],
                data[field_index['auth_seq_id']],
                x, y, z,
                data[field_index['id']],
                data[field_index['type_symbol']]
            ]
            
            if residue_name not in ligands:
                ligands[residue_name] = []
            ligands[residue_name].append(atom_fields)
        
        return ligands, None

    def process_single_cif(cif_file):
        pdb_id = os.path.basename(cif_file).split('.')[0]
        try:
            with open(cif_file, 'r', encoding='utf-8') as f:
                cif_content = f.read()
            
            ligands, err = parse_cif_atom_sites(cif_content)
            if err:
                raise ValueError(err)
            
            for residue_name, atom_fields_list in ligands.items():
                ligand_cif_file = os.path.join(ligands_dir, f'{pdb_id}_{residue_name}.cif')
                
                cif_lines = [
                    f"data_{pdb_id}_{residue_name}\n",
                    "loop_\n",
                    "_atom_site.group_PDB\n",
                    "_atom_site.auth_comp_id\n",
                    "_atom_site.auth_asym_id\n",
                    "_atom_site.auth_seq_id\n",
                    "_atom_site.Cartn_x\n",
                    "_atom_site.Cartn_y\n",
                    "_atom_site.Cartn_z\n",
                    "_atom_site.id\n",
                    "_atom_site.type_symbol\n"
                ]
                for atom_fields in atom_fields_list:
                    cif_lines.append(
                        f"{atom_fields[0]} {atom_fields[1]} {atom_fields[2]} {atom_fields[3]} "
                        f"{atom_fields[4]:.3f} {atom_fields[5]:.3f} {atom_fields[6]:.3f} "
                        f"{atom_fields[7]} {atom_fields[8]}\n"
                    )
                
                with open(ligand_cif_file, 'w', encoding='utf-8') as f:
                    f.writelines(cif_lines)
                print(f"成功保存配体 {residue_name} 的CIF文件: {ligand_cif_file}")

        except Exception as e:
            error_msg = f"处理 {cif_file} 时出错: {e}"
            print(error_msg)
            with open(error_file, 'a', encoding='utf-8') as f:
                f.write(error_msg + '\n')

    with ThreadPoolExecutor() as executor:
        executor.map(process_single_cif, cif_files)

def pdb_main(args):
    exclude_residues = parse_exclude_residues()
    
    for section in pdb_config.sections():
        output_dir = pdb_config.get(section, 'output')
        uniprot_id = pdb_config.get(section, 'uniprot')

        pdb_ids = get_pdb_ids_from_uniprot(uniprot_id)
        if pdb_ids:
            download_pdb_files(pdb_ids, output_dir)
        else:
            print(f"任务 {section}: 未找到对应的PDB ID。")

        ligand_pdb_dict = extract_ligands_from_pdb(output_dir, exclude_residues)
        write_ligand_info_to_md(ligand_pdb_dict, output_dir)
        move_no_ligand_files(output_dir, ligand_pdb_dict)
        move_with_ligand_files(output_dir, ligand_pdb_dict)

        with open(os.path.join(output_dir, 'pdb_ligand.md'), 'r', encoding='utf-8') as f:
            lines = f.readlines()
        ligands = []
        for line in lines[2:]:
            parts = line.strip().split('|')
            if len(parts) >= 2:
                ligand_str = parts[1].strip()
                ligands.append(ligand_str)
        unique_ligands = set(ligands)

        download_ligand_json(unique_ligands, output_dir)
        write_chemical_info_to_md(unique_ligands, output_dir)
        extract_ligands_coordinates(output_dir, exclude_residues)

    print('任务完成。')

# ------------------------------ 主程序入口 ------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="蛋白质分析综合工具")
    subparsers = parser.add_subparsers(dest="tool", required=True)

    # 序列分析工具子解析器
    seq_parser = subparsers.add_parser("seq", help="蛋白质序列处理工具")
    seq_subparsers = seq_parser.add_subparsers(dest="command", required=True)

    # seq fetch命令
    fetch_parser = seq_subparsers.add_parser("fetch", help="批量获取蛋白质序列和结构域信息")
    fetch_parser.add_argument("input_file", help="包含HGNC基因名称的文本文件")
    fetch_parser.add_argument("output_dir", default="protein_sequences", help="输出目录")

    # seq extract命令
    extract_parser = seq_subparsers.add_parser("extract", help="提取指定位置的序列")
    extract_parser.add_argument("fasta_file", help="输入FASTA文件路径")
    extract_parser.add_argument("start", type=int, help="起始位置")
    extract_parser.add_argument("end", type=int, help="结束位置")

    # seq mutate命令
    mut_parser = seq_subparsers.add_parser("mutate", help="执行蛋白质序列突变")
    mut_parser.add_argument("fasta_file", help="输入FASTA文件路径")
    mut_parser.add_argument("--pos", nargs="+", type=int, required=True, help="突变位置")
    mut_parser.add_argument("--aa", nargs="+", required=True, help="新氨基酸")

    # seq align命令
    align_parser = seq_subparsers.add_parser("align", help="比较多个序列")
    align_parser.add_argument("fasta_files", nargs="+", help="FASTA文件路径列表")

    # PDB分析工具子解析器
    pdb_parser = subparsers.add_parser("pdb", help="PDB文件处理工具")

    args = parser.parse_args()

    if args.tool == "seq":
        seq_main(args)
    elif args.tool == "pdb":
        pdb_main(args)