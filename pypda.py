import argparse
import os
import requests
import configparser
import json
import warnings
import shutil
from pathlib import Path
import xml.etree.ElementTree as ET
from typing import List, Dict, Tuple, Generator, Optional, Set, Any, TypedDict
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio import ExPASy
from Bio import SwissProt
from Bio.PDB.PDBList import PDBList
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from glob import glob
from concurrent.futures import ThreadPoolExecutor
from functools import lru_cache
from datetime import datetime
import sys

# 过滤 PDB 构建警告
warnings.filterwarnings("ignore", category=PDBConstructionWarning)


class ConfigManager:
    """配置文件管理类，负责读取和管理配置信息"""
    def __init__(self):
        self.seq_config = self._read_config('seq_config.ini')
        self.pdb_config = self._read_config('pdb_config.ini')
        self.ns = dict(self.seq_config.items('XML_NAMESPACES')) if self.seq_config.has_section('XML_NAMESPACES') else {}
        self.error_log_path = os.path.join(os.getcwd(), 'error.txt')
        self.uniprot_api = self.seq_config.get('UNIPROT_API', 'UNIPROT_API', fallback='https://www.uniprot.org/uniprot/')
        self.uniprot_xml_api = self.seq_config.get('UNIPROT_API', 'UNIPROT_XML_API', fallback='https://www.uniprot.org/uniprot/')
        self.uniprot_api_base_url = self.seq_config.get('uniprot', 'uniprot_api_base_url', fallback='https://rest.uniprot.org/uniprotkb/')

    @staticmethod
    def _read_config(config_path: str) -> configparser.ConfigParser:
        """读取配置文件

        Args:
            config_path: 配置文件路径

        Returns:
            配置解析器对象
        """
        config = configparser.ConfigParser()
        config.read(config_path)
        return config


class Logger:
    """日志管理类，负责统一错误处理和日志输出"""
    @staticmethod
    def log_error(error_msg: str, error_file_path: str) -> None:
        """记录错误信息到日志文件并打印

        Args:
            error_msg: 错误信息
            error_file_path: 错误日志文件路径
        """
        with open(error_file_path, 'a', encoding='utf-8') as f:
            f.write(f"{error_msg}\n")
        print(error_msg)


class CommonUtils:
    """通用工具类，提供项目中常用的工具函数"""
    @staticmethod
    def read_fasta(file_path: str, return_header: bool = False) -> Tuple[str, str] or str:
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

    @staticmethod
    def generate_output_filename(base_name: str, ext: str, *args: Any) -> str:
        """统一生成输出文件名

        Args:
            base_name: 基础文件名
            ext: 文件扩展名（含.）
            *args: 可变参数，用于生成附加信息

        Returns:
            生成的文件名
        """
        if not args:
            return f"{base_name}{ext}"
        info = "_".join(map(str, args)) if len(args) > 1 else str(args[0])
        return f"{base_name}_{info}{ext}"

    @staticmethod
    def save_sequence(header: str, sequence: str, output_file: str) -> None:
        """保存序列到FASTA文件

        Args:
            header: 序列头部信息
            sequence: 序列内容
            output_file: 输出文件路径
        """
        with open(output_file, 'w', encoding="utf-8") as file:
            file.write(header + '\n')
            for i in range(0, len(sequence), 60):
                file.write(sequence[i:i+60] + '\n')

    @staticmethod
    def parallel_executor(func: callable, items: List[Any], max_workers: Optional[int] = None) -> None:
        """并行执行函数

        Args:
            func: 要并行执行的函数
            items: 迭代参数列表
            max_workers: 最大工作线程数
        """
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            executor.map(func, items)

    @staticmethod
    def save_to_json(data: Dict[str, Any], accession: str) -> str:
        """将数据保存为JSON文件

        Args:
            data: 要保存的字典数据
            accession: UniProt编号，用于生成文件名

        Returns:
            保存的JSON文件路径
        """
        today = datetime.now().strftime("%Y%m%d")
        json_filename = f'{accession}_{today}.json'
        try:
            with open(json_filename, 'w', encoding='utf-8') as f:
                json.dump(data, f, ensure_ascii=False, indent=4)
            print(f"JSON数据已保存至 {json_filename}")
            return json_filename
        except IOError as e:
            Logger.log_error(f"保存JSON文件失败: {str(e)}", Logger.error_log_path)
            raise

    @staticmethod
    def load_from_json(file_path: str) -> Dict[str, Any]:
        """从JSON文件加载数据

        Args:
            file_path: JSON文件路径

        Returns:
            从文件加载的字典数据
        """
        if not os.path.exists(file_path):
            Logger.log_error(f"文件 {file_path} 不存在", Logger.error_log_path)
            raise FileNotFoundError(f"文件 {file_path} 不存在")

        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                return json.load(f)
        except json.JSONDecodeError as e:
            Logger.log_error(f"JSON文件解析错误: {str(e)}", Logger.error_log_path)
            raise
        except IOError as e:
            Logger.log_error(f"读取文件失败: {str(e)}", Logger.error_log_path)
            raise


class SequenceProcessor:
    """序列处理类，负责蛋白质序列相关操作"""
    def __init__(self, config: ConfigManager, logger: Logger):
        self.config = config
        self.logger = logger

    def fetch_protein_sequences(self, genes: List[str], output_dir: str) -> Generator[Tuple[str, str], None, None]:
        """从UniProt获取蛋白质序列

        Args:
            genes: 基因名称列表
            output_dir: 输出目录

        Yields:
            基因名称和对应的UniProt ID元组
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        for gene in genes:
            params = {
                "query": f"gene_exact:{gene} AND organism_id:9606",
                "format": "fasta",
                "fields": "sequence"
            }
            try:
                response = requests.get(self.config.uniprot_api, params=params)
                response.raise_for_status()
                if not response.text:
                    self.logger.log_error(f"未找到基因 {gene} 对应的蛋白质序列", self.config.error_log_path)
                    continue

                output_file = output_dir / f"{gene}.fasta"
                with open(output_file, "w", encoding="utf-8") as f:
                    first_seq = response.text.split('>', 2)[1] if '>' in response.text else response.text
                    f.write(f">{first_seq}")

                print(f"已将 {gene} 的第一个全长蛋白质序列保存到 {output_file}")
                uniprot_id = first_seq.split('\n', 1)[0].split('|')[1]
                yield gene, uniprot_id

            except requests.exceptions.RequestException as e:
                self.logger.log_error(f"获取 {gene} 失败: {str(e)}", self.config.error_log_path)

    def fetch_domain_information(self, gene_uniprot_pairs: List[Tuple[str, str]], output_dir: str) -> None:
        """获取蛋白质结构域信息并保存

        Args:
            gene_uniprot_pairs: 基因和UniProt ID的元组列表
            output_dir: 输出目录
        """
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

                xml_response = requests.get(f"{self.config.uniprot_xml_api}{uniprot_id}.xml")
                xml_response.raise_for_status()
                root = ET.fromstring(xml_response.text)

                with open(domain_info_file, "a", encoding="utf-8") as domain_f:
                    domain_f.write(f"## {gene}\n")
                    domain_f.write(f"### 总氨基酸数量: {total_amino_acids}\n\n")

                    for feature in root.findall('.//uniprot:feature[@type="domain"]', self.config.ns):
                        begin = int(feature.find('uniprot:location/uniprot:begin', self.config.ns).attrib['position'])
                        end = int(feature.find('uniprot:location/uniprot:end', self.config.ns).attrib['position'])
                        domain_name = feature.attrib.get('description', '未知结构域')
                        domain_f.write(f"- 结构域名称: {domain_name}, 序列范围: {begin}-{end}\n")
                        print(f"{gene} 的 {domain_name} 结构域的序列编号范围: {begin}-{end}")

            except requests.exceptions.RequestException as e:
                self.logger.log_error(f"获取 {gene} 的结构域信息失败: {str(e)}", self.config.error_log_path)

    @staticmethod
    def extract_subsequence(fasta_file: str, start: int, end: int) -> Optional[str]:
        """提取序列的指定区域

        Args:
            fasta_file: FASTA文件路径
            start: 起始位置
            end: 结束位置

        Returns:
            提取的子序列或None（如果提取位置无效）
        """
        header, sequence = CommonUtils.read_fasta(fasta_file, return_header=True)
        if start < 1 or end > len(sequence):
            print('提取位置超出序列范围')
            return None
        return sequence[start - 1:end]

    @staticmethod
    def validate_input(mutation_positions: List[int], new_amino_acids: List[str]) -> bool:
        """验证突变输入的有效性

        Args:
            mutation_positions: 突变位置列表
            new_amino_acids: 新氨基酸列表

        Returns:
            如果输入有效则返回True，否则退出程序
        """
        if len(mutation_positions) != len(new_amino_acids):
            print("突变位置和新氨基酸的数量必须相同。")
            exit(1)
        return True

    @staticmethod
    def perform_mutations(record: Any, mutation_positions: List[int], new_amino_acids: List[str]) -> Any:
        """执行序列突变

        Args:
            record: SeqRecord对象
            mutation_positions: 突变位置列表
            new_amino_acids: 新氨基酸列表

        Returns:
            突变后的SeqRecord对象
        """
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord

        sequence = str(record.seq)
        for pos, aa in zip(mutation_positions, new_amino_acids):
            sequence = sequence[:pos - 1] + aa + sequence[pos:]
        return SeqRecord(Seq(sequence), id=record.id, description=record.description)

    @staticmethod
    def read_align_fasta(file_path: str) -> str:
        """读取比对用的FASTA文件

        Args:
            file_path: 文件路径

        Returns:
            序列字符串
        """
        return str(SeqIO.read(file_path, "fasta").seq)

    @staticmethod
    def pairwise_alignment(seq1: str, seq2: str) -> Tuple[float, float]:
        """执行双序列比对并计算同源性

        Args:
            seq1: 第一条序列
            seq2: 第二条序列

        Returns:
            比对得分和同源性百分比
        """
        aligner = PairwiseAligner()
        best_alignment = next(aligner.align(seq1, seq2))
        score = best_alignment.score
        homology = (score / max(len(seq1), len(seq2))) * 100
        return score, homology

    def compare_sequences(self, file_paths: List[str]) -> None:
        """比较多个FASTA文件中的序列

        Args:
            file_paths: FASTA文件路径列表
        """
        sequences = [self.read_align_fasta(fp) for fp in file_paths]
        for i, seq1 in enumerate(sequences):
            for j, seq2 in enumerate(sequences[i+1:], i+1):
                score, homology = self.pairwise_alignment(seq1, seq2)
                print(f"比对 {i+1} 和 {j+1}:\n比对得分: {score}\n同源性百分比: {homology:.2f}%\n")

    def process_command(self, args: argparse.Namespace) -> None:
        """处理序列相关命令

        Args:
            args: 命令行参数
        """
        if args.command == "fetch":
            with open(args.input_file) as f:
                genes = [line.strip() for line in f if line.strip()]
            if genes:
                gene_uniprot_pairs = list(self.fetch_protein_sequences(genes, args.output_dir))
                self.fetch_domain_information(gene_uniprot_pairs, args.output_dir)
            else:
                print("输入文件中未找到有效的基因名称")

        elif args.command == "extract":
            if not os.path.isfile(args.fasta_file):
                print(f"错误：文件 {args.fasta_file} 不存在")
                return
            header, _ = CommonUtils.read_fasta(args.fasta_file, return_header=True)
            extracted_seq = self.extract_subsequence(args.fasta_file, args.start, args.end)
            output_file = CommonUtils.generate_output_filename(args.fasta_file, ".fasta", args.start, args.end)
            CommonUtils.save_sequence(header, extracted_seq, output_file)
            print(f"已保存截取序列至 {output_file}")

        elif args.command == "mutate":
            self.validate_input(args.pos, args.aa)
            record = next(SeqIO.parse(args.fasta_file, "fasta"))
            mutated_record = self.perform_mutations(record, args.pos, args.aa)
            mutation_info = "".join([f"{pos}{aa}" for pos, aa in zip(args.pos, args.aa)])
            output_file = CommonUtils.generate_output_filename(os.path.splitext(args.fasta_file)[0], ".fasta", mutation_info)
            SeqIO.write(mutated_record, output_file, "fasta")
            print(f"已保存突变序列至 {output_file}")

        elif args.command == "align":
            self.compare_sequences(args.fasta_files)


class PDBProcessor:
    """PDB处理类，负责PDB文件相关操作"""
    def __init__(self, config: ConfigManager, logger: Logger):
        self.config = config
        self.logger = logger
        self.exclude_residues = self.parse_exclude_residues()

    def parse_exclude_residues(self) -> List[str]:
        """解析排除残基的配置文件

        Returns:
            需要排除的残基列表
        """
        config = configparser.ConfigParser()
        config.read('exclude_residues.ini')
        exclude_residues = []
        for section in config.sections():
            for key in config[section]:
                exclude_residues.extend([res.strip() for res in config[section][key].split(',')])
        return exclude_residues

    @lru_cache(maxsize=None)
    def get_pdb_ids_from_uniprot(self, uniprot_id: str) -> List[str]:
        """通过UniProt ID获取对应的PDB ID列表（带缓存优化）

        Args:
            uniprot_id: UniProt ID

        Returns:
            PDB ID列表
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
            self.logger.log_error(f"获取PDB ID时出错: {e}", self.config.error_log_path)
            return []

    def download_pdb_files(self, pdb_ids: List[str], output_dir: str) -> None:
        """下载PDB文件

        Args:
            pdb_ids: PDB ID列表
            output_dir: 输出目录
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        pdbl = PDBList()
        failed_files = []
        error_file = output_dir / 'error.txt'

        def download_single_file(pdb_id: str) -> None:
            try:
                pdbl.retrieve_pdb_file(pdb_id, pdir=output_dir, file_format='mmCif')
                print(f"成功下载CIF文件: {pdb_id}.cif")
            except Exception as e:
                error_msg = f"下载CIF文件 {pdb_id} 时出错: {e}"
                print(error_msg)
                failed_files.append(pdb_id)
                with open(error_file, 'a', encoding='utf-8') as f:
                    f.write(error_msg + '\n')

        CommonUtils.parallel_executor(download_single_file, pdb_ids)

        if failed_files:
            with open(os.path.join(output_dir, 'download_failed.txt'), 'a', encoding='utf-8') as f:
                for failed_file in failed_files:
                    f.write(f"{failed_file}.cif\n")

    def extract_ligands_from_pdb(self, folder_path: str) -> Dict[str, List[str]]:
        """从.cif文件中提取配体信息

        Args:
            folder_path: 文件夹路径

        Returns:
            配体到PDB ID的映射字典
        """
        parser = MMCIFParser()
        ligand_pdb_dict: Dict[str, List[str]] = {}
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
                                    if residue_name not in self.exclude_residues and residue_name not in ligands:
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

    @staticmethod
    def write_ligand_info_to_md(ligand_pdb_dict: Dict[str, List[str]], output_dir: str) -> None:
        """将配体信息写入pdb_ligand.md文件

        Args:
            ligand_pdb_dict: 配体到PDB ID的映射字典
            output_dir: 输出目录
        """
        output_file = 'pdb_ligand.md'
        with open(os.path.join(output_dir, output_file), 'w') as outfile:
            outfile.write('| Ligands | PDB ID |\n')
            outfile.write('| --- | --- |\n')
            for ligand, pdb_ids in ligand_pdb_dict.items():
                pdb_id_str = ', '.join(pdb_ids)
                outfile.write(f'| {ligand} | {pdb_id_str} |\n')
        print(f'结果已写入 {os.path.join(output_dir, output_file)}')

    def move_no_ligand_files(self, folder_path: str, ligand_pdb_dict: Dict[str, List[str]]) -> None:
        """将没有配体的CIF文件移动到no_ligand文件夹

        Args:
            folder_path: 文件夹路径
            ligand_pdb_dict: 配体到PDB ID的映射字典
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
                        self.logger.log_error(f"移动 {filename} 时出错: {e}", error_file)

    def move_with_ligand_files(self, folder_path: str, ligand_pdb_dict: Dict[str, List[str]]) -> None:
        """将含有配体的CIF文件移动到with_ligands文件夹

        Args:
            folder_path: 文件夹路径
            ligand_pdb_dict: 配体到PDB ID的映射字典
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
                        self.logger.log_error(f"移动含配体文件 {filename} 时出错: {e}", error_file)

    def download_ligand_json(self, unique_ligands: Set[str], output_dir: str) -> None:
        """通过API并行查询配体的JSON文件并下载保存到json子文件夹

        Args:
            unique_ligands: 唯一配体集合
            output_dir: 输出目录
        """
        json_dir = os.path.join(output_dir, 'json')
        os.makedirs(json_dir, exist_ok=True)

        def download_single_ligand(ligand: str) -> None:
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

        CommonUtils.parallel_executor(download_single_ligand, list(unique_ligands))

    def write_chemical_info_to_md(self, unique_ligands: Set[str], output_dir: str) -> None:
        """从json子文件夹中读取JSON文件并写入chemical_components_info.md

        Args:
            unique_ligands: 唯一配体集合
            output_dir: 输出目录
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

    def extract_ligands_coordinates(self, output_dir: str) -> None:
        """直接通过文本解析.cif文件提取配体坐标并保存成cif文件（并行版）

        Args:
            output_dir: 输出目录
        """
        with_ligands_dir = os.path.join(output_dir, 'with_ligands')
        ligands_dir = os.path.join(with_ligands_dir, 'ligands')
        os.makedirs(ligands_dir, exist_ok=True)

        cif_files = glob(os.path.join(with_ligands_dir, '*.cif'))
        error_file = os.path.join(output_dir, 'error.txt')

        def parse_cif_atom_sites(cif_content: str) -> Tuple[Dict[str, List[List[Any]]], Optional[str]]:
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
                return {}, "缺少必要的_atom_site字段"

            ligands = {}
            for data in atom_site_data:
                if data[field_index['group_PDB']] != 'HETATM':
                    continue

                residue_name = data[field_index['auth_comp_id']].strip()
                if residue_name in self.exclude_residues:
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

        def process_single_cif(cif_file: str) -> None:
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

        CommonUtils.parallel_executor(process_single_cif, cif_files)

    def process(self) -> None:
        """处理PDB相关任务"""
        for section in self.config.pdb_config.sections():
            output_dir = self.config.pdb_config.get(section, 'output')
            uniprot_id = self.config.pdb_config.get(section, 'uniprot')

            pdb_ids = self.get_pdb_ids_from_uniprot(uniprot_id)
            if pdb_ids:
                self.download_pdb_files(pdb_ids, output_dir)
            else:
                print(f"任务 {section}: 未找到对应的PDB ID。")

            ligand_pdb_dict = self.extract_ligands_from_pdb(output_dir)
            self.write_ligand_info_to_md(ligand_pdb_dict, output_dir)
            self.move_no_ligand_files(output_dir, ligand_pdb_dict)
            self.move_with_ligand_files(output_dir, ligand_pdb_dict)

            with open(os.path.join(output_dir, 'pdb_ligand.md'), 'r', encoding='utf-8') as f:
                lines = f.readlines()
            ligands = []
            for line in lines[2:]:
                parts = line.strip().split('|')
                if len(parts) >= 2:
                    ligand_str = parts[1].strip()
                    ligands.append(ligand_str)
            unique_ligands = set(ligands)

            self.download_ligand_json(unique_ligands, output_dir)
            self.write_chemical_info_to_md(unique_ligands, output_dir)
            self.extract_ligands_coordinates(output_dir)

        print('任务完成。')


class ProteinInfo(TypedDict):
    """蛋白质信息结构化数据类型定义

    包含从UniProt JSON数据中提取的各类蛋白质信息
    """
    basic_info: Dict[str, Any]
    biology_info: Dict[str, Any]
    protein_desc: Dict[str, Any]
    comments: Dict[str, List[Dict[str, Any]]]
    features: Dict[str, Any]
    interactions: List[Dict[str, Any]]
    keywords: List[str]
    references: List[Dict[str, Any]]
    cross_references: Dict[str, List[Dict[str, Any]]]
    sequence: Dict[str, Any]


class UniProtAPI:
    """UniProt API交互类，负责从UniProt API获取数据"""

    def __init__(self, config: ConfigManager, logger: Logger):
        """初始化UniProtAPI实例

        Args:
            config: 配置管理器实例
            logger: 日志记录器实例
        """
        self.config = config
        self.logger = logger
        self.api_base_url = config.uniprot_api_base_url
        self.session = requests.Session()
        self.session.headers.update({
            'Accept': 'application/json',
            'User-Agent': 'UniProtAnalyzer/1.0'
        })

    def get_uniprot_data(self, accession: str) -> Optional[Dict[str, Any]]:
        """根据UniProt编号从API获取数据"""
        url = f'{self.api_base_url}/{accession}'
        try:
            self.logger.log_error(f"正在从 {url} 获取数据...", self.config.error_log_path)
            response = self.session.get(url, timeout=10)
            response.raise_for_status()
            return response.json()
        except requests.exceptions.HTTPError as e:
            self.logger.log_error(f"HTTP请求错误: {str(e)}", self.config.error_log_path)
        except requests.exceptions.ConnectionError:
            self.logger.log_error("网络连接错误", self.config.error_log_path)
        except requests.exceptions.Timeout:
            self.logger.log_error("请求超时", self.config.error_log_path)
        except Exception as e:
            self.logger.log_error(f"获取UniProt数据失败: {str(e)}", self.config.error_log_path)
        return None


class ProteinAnalyzer:
    """蛋白质数据分析类，负责从JSON数据中提取和组织蛋白质信息"""

    @staticmethod
    def extract_protein_info(data: Dict[str, Any]) -> ProteinInfo:
        """从JSON数据中提取蛋白质信息并组织为结构化字典

        Args:
            data: 从UniProt API获取的原始JSON数据

        Returns:
            包含各类蛋白质信息的结构化字典
        """
        info: ProteinInfo = {
            'basic_info': {},
            'biology_info': {},
            'protein_desc': {},
            'comments': {},
            'features': {},
            'interactions': [],
            'keywords': [],
            'references': [],
            'cross_references': {},
            'sequence': {}
        }

        # 1. 基本识别信息
        info['basic_info'] = {
            'entryType': data.get('entryType', 'N/A'),
            'primaryAccession': data.get('primaryAccession', 'N/A'),
            'secondaryAccessions': data.get('secondaryAccessions', []),
            'uniProtkbId': data.get('uniProtkbId', 'N/A'),
            'annotationScore': data.get('annotationScore', 'N/A'),
            'entryAudit': data.get('entryAudit', {})
        }

        # 2. 生物学背景信息
        organism = data.get('organism', {})
        info['biology_info'] = {
            'scientificName': organism.get('scientificName', 'N/A'),
            'commonName': organism.get('commonName', 'N/A'),
            'taxonId': organism.get('taxonId', 'N/A'),
            'lineage': organism.get('lineage', []),
            'proteinExistence': data.get('proteinExistence', 'N/A')
        }

        # 3. 蛋白质描述和基因信息
        protein_desc = data.get('proteinDescription', {})
        recommended_name = protein_desc.get('recommendedName', {})
        info['protein_desc'] = {
            'recommendedName': {
                'fullName': recommended_name.get('fullName', {}).get('value', 'N/A'),
                'shortNames': [sn.get('value') for sn in recommended_name.get('shortNames', []) if sn.get('value')],
                'ecNumbers': [ec.get('value') for ec in recommended_name.get('ecNumbers', []) if ec.get('value')]
            },
            'alternativeNames': protein_desc.get('alternativeNames', []),
            'genes': data.get('genes', [])
        }

        # 4. 功能和活性注释 - 修改为大小写不敏感匹配
        comments = data.get('comments', [])
        comment_types = ['FUNCTION', 'CATALYTIC ACTIVITY', 'COFACTOR', 'ACTIVITY REGULATION',
                         'TISSUE SPECIFICITY', 'SUBCELLULAR LOCATION', 'PTM', 'SIMILARITY']
        # 创建大小写不敏感的映射字典
        comment_type_map = {ct.lower(): ct for ct in comment_types}
        info['comments'] = {ct: [] for ct in comment_types}
        
        for comment in comments:
            # 将API返回的commentType转换为小写进行匹配
            ct_lower = comment.get('commentType', '').lower()
            if ct_lower in comment_type_map:
                ct = comment_type_map[ct_lower]
                info['comments'][ct].append(comment)

        # 5. 蛋白质特征
        features = data.get('features', [])
        info['features'] = {
            'counts': {},
            'detailed': {}
        }
        detailed_feature_types = ["Chain", "Region", "Active site", "Binding site", "Modified residue", "Mutagenesis"]
        for dt in detailed_feature_types:
            info['features']['detailed'][dt] = []

        for feature in features:
            ft = feature.get('type')
            info['features']['counts'][ft] = info['features']['counts'].get(ft, 0) + 1
            if ft in detailed_feature_types:
                info['features']['detailed'][ft].append(feature)

        # 6. 蛋白质相互作用
        info['interactions'] = data.get('interactions', [])

        # 7. 关键词
        info['keywords'] = [kw.get('name') for kw in data.get('keywords', []) if kw.get('name')]

        # 8. 参考文献
        info['references'] = data.get('references', [])

        # 9. 交叉引用
        cross_references = data.get('uniProtKBCrossReferences', [])
        important_databases = ["PDB", "DrugBank", "GO", "Reactome", "HGNC", "GeneID", "KEGG"]
        info['cross_references'] = {db: [] for db in important_databases}
        for xref in cross_references:
            db = xref.get('database')
            if db in info['cross_references']:
                info['cross_references'][db].append(xref)

        # 10. 蛋白质序列信息
        info['sequence'] = data.get('sequence', {})

        return info


class ReportGenerator:
    """报告生成器类，负责将蛋白质信息生成为Markdown报告"""

    @staticmethod
    def generate_md_report(info: ProteinInfo, md_filename: str) -> None:
        """将分析结果生成Markdown文件

        Args:
            info: 包含蛋白质信息的结构化字典
            md_filename: 要生成的Markdown文件名
        """
        try:
            with open(md_filename, 'w', encoding='utf-8') as f:
                f.write("# UniProt 蛋白质信息分析报告\n\n")

                # 1. 基本识别信息
                f.write("## 1. 基本识别信息\n")
                bi = info['basic_info']
                f.write(f"- **条目类型**: {bi['entryType']}\n")
                f.write(f"- **主要登录号**: {bi['primaryAccession']}\n")
                if bi['secondaryAccessions']:
                    f.write(f"- **次要登录号**: {', '.join(bi['secondaryAccessions'])}\n")
                f.write(f"- **UniProtKB ID**: {bi['uniProtkbId']}\n")
                f.write(f"- **注释评分**: {bi['annotationScore']}\n")
                audit = bi['entryAudit']
                f.write(f"- **首次公开日期**: {audit.get('firstPublicDate', 'N/A')}\n")
                f.write(f"- **最后注释更新日期**: {audit.get('lastAnnotationUpdateDate', 'N/A')}\n")
                f.write(f"- **最后序列更新日期**: {audit.get('lastSequenceUpdateDate', 'N/A')}\n")
                f.write(f"- **条目版本**: {audit.get('entryVersion', 'N/A')}\n")
                f.write(f"- **序列版本**: {audit.get('sequenceVersion', 'N/A')}\n\n")

                # 2. 生物学背景信息
                f.write("## 2. 生物学背景信息\n")
                bio = info['biology_info']
                f.write(f"- **科学名称**: {bio['scientificName']}\n")
                f.write(f"- **常用名称**: {bio['commonName']}\n")
                f.write(f"- **分类ID**: {bio['taxonId']}\n")
                if bio['lineage']:
                    f.write(f"- **生物学谱系**: {' -> '.join(bio['lineage'])}\n")
                f.write(f"- **蛋白质存在证据**: {bio['proteinExistence']}\n\n")

                # 3. 蛋白质描述和基因信息
                f.write("## 3. 蛋白质描述和基因信息\n")
                pd = info['protein_desc']
                rn = pd['recommendedName']
                f.write(f"- **推荐全名**: {rn['fullName']}\n")
                if rn['shortNames']:
                    f.write(f"- **推荐简称**: {', '.join(rn['shortNames'])}\n")
                if rn['ecNumbers']:
                    f.write(f"- **EC 编号**: {', '.join(rn['ecNumbers'])}\n")

                alternative_names = pd['alternativeNames']
                if alternative_names:
                    f.write("- **备选名称**:\n")
                    for alt_name in alternative_names:
                        alt_full_name = alt_name.get('fullName', {}).get('value', 'N/A')
                        alt_ec_numbers = [ec.get('value') for ec in alt_name.get('ecNumbers', []) if ec.get('value')]
                        f.write(f"  - {alt_full_name}" + (f" (EC: {', '.join(alt_ec_numbers)})" if alt_ec_numbers else "") + "\n")

                genes = pd['genes']
                if genes:
                    for gene in genes:
                        gene_name = gene.get('geneName', {}).get('value', 'N/A')
                        f.write(f"- **基因名称**: {gene_name}\n")
                        gene_synonyms = [syn.get('value') for syn in gene.get('synonyms', []) if syn.get('value')]
                        if gene_synonyms:
                            f.write(f"  - **基因同义词**: {', '.join(gene_synonyms)}\n")
                f.write("\n")

                # 4. 功能和活性注释 - 增强文本提取逻辑
                f.write("## 4. 功能和活性注释\n")
                comments = info['comments']
                # 确保按指定顺序输出注释类型
                for ct in ['FUNCTION', 'CATALYTIC ACTIVITY', 'COFACTOR', 'ACTIVITY REGULATION',
                           'TISSUE SPECIFICITY', 'SUBCELLULAR LOCATION', 'PTM', 'SIMILARITY']:
                    if comments[ct]:
                        f.write(f"### {ct}\n")
                        for comment in comments[ct]:
                            # 处理可能的嵌套text结构
                            texts = []
                            for text in comment.get('texts', []):
                                if isinstance(text, dict) and 'value' in text:
                                    texts.append(text['value'])
                                elif isinstance(text, str):
                                    texts.append(text)
                            if texts:
                                f.write(f"- {'; '.join(texts)}\n")
                f.write("\n")

                # 5. 蛋白质特征
                f.write("## 5. 蛋白质特征\n")
                features = info['features']
                f.write("### 特征类型统计\n")
                for ft, count in features['counts'].items():
                    f.write(f"- {ft}: {count}\n")

                f.write("\n### 详细特征信息\n")
                for dt in features['detailed']:
                    if features['detailed'][dt]:
                        f.write(f"#### {dt}\n")
                        for feature in features['detailed'][dt]:
                            desc = feature.get('description', 'N/A')
                            if isinstance(desc, dict):
                                desc = desc.get('value', 'N/A')

                            loc = feature.get('location', {})
                            begin = loc.get('start', {}).get('value', 'N/A')
                            end = loc.get('end', {}).get('value', 'N/A')
                            f.write(f"- {desc} (位置: {begin}-{end})\n")
                f.write("\n")

                # 6. 蛋白质相互作用
                f.write("## 6. 蛋白质相互作用\n")
                interactions = info['interactions']
                if interactions:
                    for i, interaction in enumerate(interactions, 1):
                        interactor = interaction.get('interactor', {})
                        interactor_id = interactor.get('uniprotId', 'N/A')
                        interactor_name = interactor.get('name', 'N/A')
                        f.write(f"- 相互作用蛋白 {i}: {interactor_name} (UniProt ID: {interactor_id})\n")
                else:
                    f.write("- 未找到相互作用信息\n")
                f.write("\n")

                # 7. 关键词
                f.write("## 7. 关键词\n")
                keywords = info['keywords']
                if keywords:
                    f.write(f"- {', '.join(keywords)}\n")
                else:
                    f.write("- 未找到关键词信息\n")
                f.write("\n")

                # 8. 参考文献
                f.write("## 8. 参考文献\n")
                references = info['references']
                if references:
                    for i, ref in enumerate(references[:5], 1):
                        citation = ref.get('citation', {})
                        authors = citation.get('authors', 'N/A')
                        title = citation.get('title', 'N/A')
                        journal = citation.get('journal', 'N/A')
                        year = citation.get('publicationDate', 'N/A')[:4]  # 提取年份
                        f.write(f"- [{i}] {authors}, {title}, {journal}, {year}\n")
                    if len(references) > 5:
                        f.write(f"- 显示前5篇，共{len(references)}篇参考文献\n")
                else:
                    f.write("- 未找到参考文献信息\n")
                f.write("\n")

                # 9. 交叉引用
                f.write("## 9. 交叉引用\n")
                cross_refs = info['cross_references']
                for db in cross_refs:
                    if cross_refs[db]:
                        f.write(f"### {db}\n")
                        ids = [xref.get('id') for xref in cross_refs[db] if xref.get('id')]
                        f.write(f"- {', '.join(ids)}\n")
                        if len(cross_refs[db]) > 5:
                            f.write(f"- 共{len(cross_refs[db])}个条目\n")
                f.write("\n")

                # 10. 蛋白质序列信息
                f.write("## 10. 蛋白质序列信息\n")
                sequence = info['sequence']
                f.write(f"- **序列长度**: {sequence.get('length', 'N/A')} 个氨基酸\n")
                f.write(f"- **序列版本**: {sequence.get('sequenceVersion', 'N/A')}\n")
                f.write(f"- **序列MD5**: {sequence.get('md5', 'N/A')}\n")
                seq = sequence.get('value', '')
                if seq:
                    # 每80个字符换行显示
                    formatted_seq = '\n'.join([seq[i:i+80] for i in range(0, len(seq), 80)])
                    f.write("- **氨基酸序列**:\n\n```\n{}\n```\n".format(formatted_seq))

            print(f"Markdown报告已生成至 {md_filename}")
        except IOError as e:
            Logger.log_error(f"生成Markdown报告失败: {str(e)}", Logger.error_log_path)
            raise


class PypdaApp:
    """Pypda应用主类"""
    def __init__(self):
        self.config = ConfigManager()
        self.logger = Logger()
        self.seq_processor = SequenceProcessor(self.config, self.logger)
        self.pdb_processor = PDBProcessor(self.config, self.logger)
        self.uniprot_api = UniProtAPI(self.config, self.logger)

    def setup_parser(self) -> argparse.ArgumentParser:
        """设置命令行参数解析器

        Returns:
            配置好的参数解析器
        """
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

        # UniProt分析工具子解析器
        uniprot_parser = subparsers.add_parser("uniprot", help="UniProt数据处理工具")
        uniprot_subparsers = uniprot_parser.add_subparsers(dest="command", required=True)

        # uniprot fetch命令
        uniprot_fetch_parser = uniprot_subparsers.add_parser("fetch", help="从UniProt API获取数据并生成报告")
        uniprot_fetch_parser.add_argument("accession", help="UniProt蛋白质编号，例如: Q13547")

        # uniprot analyze命令
        uniprot_analyze_parser = uniprot_subparsers.add_parser("analyze", help="分析现有JSON文件并生成报告")
        uniprot_analyze_parser.add_argument("-f", "--file", required=True, help="要分析的JSON文件路径")

        return parser

    def run(self) -> None:
        """运行应用程序"""
        parser = self.setup_parser()
        args = parser.parse_args()

        if args.tool == "seq":
            self.seq_processor.process_command(args)
        elif args.tool == "pdb":
            self.pdb_processor.process()
        elif args.tool == "uniprot":
            if args.command == "fetch":
                # 执行fetch命令
                data = self.uniprot_api.get_uniprot_data(args.accession)
                if not data:
                    self.logger.log_error("无法获取UniProt数据", self.config.error_log_path)
                    return
                json_filename = CommonUtils.save_to_json(data, args.accession)
                analyzer = ProteinAnalyzer()
                protein_info = analyzer.extract_protein_info(data)
                md_filename = os.path.splitext(json_filename)[0] + '.md'
                ReportGenerator.generate_md_report(protein_info, md_filename)
            elif args.command == "analyze":
                # 执行analyze命令
                try:
                    data = CommonUtils.load_from_json(args.file)
                    analyzer = ProteinAnalyzer()
                    protein_info = analyzer.extract_protein_info(data)
                    md_filename = os.path.splitext(args.file)[0] + '.md'
                    ReportGenerator.generate_md_report(protein_info, md_filename)
                except Exception as e:
                    self.logger.log_error(f"分析文件失败: {str(e)}", self.config.error_log_path)


if __name__ == "__main__":
    app = PypdaApp()
    app.run()