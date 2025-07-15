import argparse
import os
import requests
import configparser
import json
import warnings
import shutil
import sys
import io
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import List, Dict, Tuple, Generator, Optional, Set, Any, TypedDict, Union
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
            f.write(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - ERROR: {error_msg}\n")
        print(f"ERROR: {error_msg}")


class CommonUtils:
    """通用工具类，提供项目中常用的工具函数"""
    @staticmethod
    def read_fasta(file_path: str, return_header: bool = False) -> Union[str, Tuple[str, str]]:
        """读取FASTA文件，可选择返回头部信息

        Args:
            file_path: 文件路径
            return_header: 是否返回头部信息

        Returns:
            序列字符串，或(头部, 序列)元组
        """
        try:
            record = SeqIO.read(file_path, "fasta")
            if return_header:
                return record.description, str(record.seq)
            return str(record.seq)
        except Exception as e:
            Logger.log_error(f"读取FASTA文件失败: {file_path} - {str(e)}", ConfigManager().error_log_path)
            raise

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
        return f"{Path(base_name).stem}_{info}{ext}" # Use Path.stem to get base name without extension

    @staticmethod
    def save_sequence(header: str, sequence: str, output_file: str) -> None:
        """保存序列到FASTA文件

        Args:
            header: 序列头部信息
            sequence: 序列内容
            output_file: 输出文件路径
        """
        try:
            with open(output_file, 'w', encoding="utf-8") as file:
                file.write(header + '\n')
                for i in range(0, len(sequence), 60):
                    file.write(sequence[i:i+60] + '\n')
            print(f"序列已保存至 {output_file}")
        except IOError as e:
            Logger.log_error(f"保存序列到文件失败: {output_file} - {str(e)}", ConfigManager().error_log_path)
            raise

    @staticmethod
    def parallel_executor(func: callable, items: List[Any], max_workers: Optional[int] = None) -> None:
        """并行执行函数

        Args:
            func: 要并行执行的函数
            items: 迭代参数列表
            max_workers: 最大工作线程数
        """
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            list(executor.map(func, items)) # Convert map object to list to ensure execution

    @staticmethod
    def save_to_json(data: Dict[str, Any], accession: str, output_dir: str = '.') -> str:
        """将数据保存为JSON文件

        Args:
            data: 要保存的字典数据
            accession: UniProt编号，用于生成文件名
            output_dir: 输出目录，默认为当前目录

        Returns:
            保存的JSON文件路径
        """
        today = datetime.now().strftime("%Y%m%d")
        json_filename = Path(output_dir) / f'{accession}_{today}.json'
        try:
            with open(json_filename, 'w', encoding='utf-8') as f:
                json.dump(data, f, ensure_ascii=False, indent=4)
            print(f"JSON数据已保存至 {json_filename}")
            return str(json_filename)
        except IOError as e:
            Logger.log_error(f"保存JSON文件失败: {json_filename} - {str(e)}", ConfigManager().error_log_path)
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
            Logger.log_error(f"文件 {file_path} 不存在", ConfigManager().error_log_path)
            raise FileNotFoundError(f"文件 {file_path} 不存在")

        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                return json.load(f)
        except json.JSONDecodeError as e:
            Logger.log_error(f"JSON文件解析错误: {file_path} - {str(e)}", ConfigManager().error_log_path)
            raise
        except IOError as e:
            Logger.log_error(f"读取文件失败: {file_path} - {str(e)}", ConfigManager().error_log_path)
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
        output_dir_path = Path(output_dir)
        output_dir_path.mkdir(parents=True, exist_ok=True)

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

                output_file = output_dir_path / f"{gene}.fasta"
                with open(output_file, "w", encoding="utf-8") as f:
                    # UniProt FASTA response might contain multiple sequences, take the first one
                    first_seq_record = next(SeqIO.parse(io.StringIO(response.text.strip()), "fasta"))
                    SeqIO.write(first_seq_record, f, "fasta")


                print(f"已将 {gene} 的第一个全长蛋白质序列保存到 {output_file}")
                # Extract UniProt ID from the header of the first sequence
                uniprot_id = first_seq_record.id.split('|')[1] if '|' in first_seq_record.id else first_seq_record.id
                yield gene, uniprot_id

            except requests.exceptions.RequestException as e:
                self.logger.log_error(f"获取 {gene} 蛋白质序列失败: {str(e)}", self.config.error_log_path)
            except StopIteration:
                self.logger.log_error(f"UniProt API返回空序列或无效FASTA格式，针对基因 {gene}", self.config.error_log_path)
            except Exception as e:
                self.logger.log_error(f"处理基因 {gene} 时发生未知错误: {str(e)}", self.config.error_log_path)


    def fetch_domain_information(self, gene_uniprot_pairs: List[Tuple[str, str]], output_dir: str) -> None:
        """获取蛋白质结构域信息并保存

        Args:
            gene_uniprot_pairs: 基因和UniProt ID的元组列表
            output_dir: 输出目录
        """
        output_dir_path = Path(output_dir)
        output_dir_path.mkdir(parents=True, exist_ok=True)
        domain_info_file = output_dir_path / "domain_info.md"

        with open(domain_info_file, "w", encoding="utf-8") as domain_f:
            domain_f.write("# 结构域信息\n\n")

        for gene, uniprot_id in gene_uniprot_pairs:
            try:
                fasta_file = output_dir_path / f"{gene}.fasta"
                total_amino_acids: Union[int, str]
                if fasta_file.exists():
                    record = SeqIO.read(fasta_file, "fasta")
                    total_amino_acids = len(record.seq)
                else:
                    total_amino_acids = "未知" # This case should ideally not happen if fetch_protein_sequences worked

                xml_response = requests.get(f"{self.config.uniprot_xml_api}{uniprot_id}.xml")
                xml_response.raise_for_status()
                root = ET.fromstring(xml_response.text)

                with open(domain_info_file, "a", encoding="utf-8") as domain_f:
                    domain_f.write(f"## {gene} (UniProt ID: {uniprot_id})\n")
                    domain_f.write(f"### 总氨基酸数量: {total_amino_acids}\n\n")
                    
                    found_domains = False
                    for feature in root.findall('.//uniprot:feature[@type="domain"]', self.config.ns):
                        begin_elem = feature.find('uniprot:location/uniprot:begin', self.config.ns)
                        end_elem = feature.find('uniprot:location/uniprot:end', self.config.ns)
                        
                        if begin_elem is not None and end_elem is not None:
                            begin = int(begin_elem.attrib['position'])
                            end = int(end_elem.attrib['position'])
                            domain_name = feature.attrib.get('description', '未知结构域')
                            domain_f.write(f"- 结构域名称: {domain_name}, 序列范围: {begin}-{end}\n")
                            print(f"{gene} 的 {domain_name} 结构域的序列编号范围: {begin}-{end}")
                            found_domains = True
                    if not found_domains:
                        domain_f.write("- 未找到结构域信息。\n")
                print(f"已将 {gene} 的结构域信息写入 {domain_info_file}")

            except requests.exceptions.RequestException as e:
                self.logger.log_error(f"获取 {gene} 的结构域信息失败: {str(e)}", self.config.error_log_path)
            except ET.ParseError as e:
                self.logger.log_error(f"解析 {gene} 的UniProt XML数据失败: {str(e)}", self.config.error_log_path)
            except Exception as e:
                self.logger.log_error(f"处理基因 {gene} 的结构域信息时发生未知错误: {str(e)}", self.config.error_log_path)


    @staticmethod
    def extract_subsequence(fasta_file: str, start: int, end: int) -> Optional[str]:
        """提取序列的指定区域

        Args:
            fasta_file: FASTA文件路径
            start: 起始位置 (1-based)
            end: 结束位置 (1-based)

        Returns:
            提取的子序列或None（如果提取位置无效）
        """
        try:
            header, sequence = CommonUtils.read_fasta(fasta_file, return_header=True)
            if not isinstance(header, str) or not isinstance(sequence, str):
                Logger.log_error(f"从文件 {fasta_file} 读取序列或头部信息失败。", ConfigManager().error_log_path)
                return None

            if start < 1 or end > len(sequence) or start > end:
                Logger.log_error(f"提取位置超出序列范围或无效 (Start: {start}, End: {end}, Seq Length: {len(sequence)})", ConfigManager().error_log_path)
                return None
            return sequence[start - 1:end]
        except Exception as e:
            Logger.log_error(f"提取子序列失败: {str(e)}", ConfigManager().error_log_path)
            return None

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
            Logger.log_error("突变位置和新氨基酸的数量必须相同。", ConfigManager().error_log_path)
            sys.exit(1)
        return True

    @staticmethod
    def perform_mutations(record: Any, mutation_positions: List[int], new_amino_acids: List[str]) -> Any:
        """执行序列突变

        Args:
            record: SeqRecord对象
            mutation_positions: 突变位置列表 (1-based)
            new_amino_acids: 新氨基酸列表

        Returns:
            突变后的SeqRecord对象
        """
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord

        sequence_list = list(str(record.seq)) # Convert to list for mutable operations
        for pos, aa in zip(mutation_positions, new_amino_acids):
            if 1 <= pos <= len(sequence_list):
                sequence_list[pos - 1] = aa
            else:
                Logger.log_error(f"突变位置 {pos} 超出序列范围。跳过此突变。", ConfigManager().error_log_path)
                # Optionally, you might want to raise an error or exit here depending on desired behavior
        mutated_sequence = "".join(sequence_list)
        return SeqRecord(Seq(mutated_sequence), id=record.id, description=record.description)

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
        # Default parameters are often fine, but can be tuned:
        # aligner.match_score = 1
        # aligner.mismatch_score = -1
        # aligner.open_gap_score = -0.5
        # aligner.extend_gap_score = -0.1

        alignments = aligner.align(seq1, seq2)
        if not alignments:
            Logger.log_error("无法找到序列比对。", ConfigManager().error_log_path)
            return 0.0, 0.0

        best_alignment = next(alignments) # Get the first (and usually best) alignment
        score = best_alignment.score
        
        # Homology calculation: identity / min(len1, len2) * 100 is often used for sequence identity
        # For overall "homology" based on alignment score, score / max_possible_score is more appropriate.
        # Max possible score for two sequences of length L1, L2 is min(L1, L2) * match_score (if aligner.match_score is 1)
        # Or, sum of individual scores of aligned residues + gap penalties.
        # A simple approach for homology based on score:
        max_len = max(len(seq1), len(seq2))
        homology = (score / max_len) * 100 if max_len > 0 else 0.0
        
        # Alternatively, for sequence identity (number of exact matches / length of shorter sequence)
        # This requires parsing the alignment itself, which PairwiseAligner doesn't directly expose as a simple count.
        # Let's stick to the current score-based homology, or refine if needed.
        
        return score, homology

    def compare_sequences(self, file_paths: List[str]) -> None:
        """比较多个FASTA文件中的序列

        Args:
            file_paths: FASTA文件路径列表
        """
        if len(file_paths) < 2:
            self.logger.log_error("至少需要两个FASTA文件才能进行比对。", self.config.error_log_path)
            return

        sequences: List[str] = []
        for fp in file_paths:
            try:
                sequences.append(self.read_align_fasta(fp))
            except Exception as e:
                self.logger.log_error(f"读取比对文件 {fp} 失败: {str(e)}", self.config.error_log_path)
                return # Stop if any file read fails

        for i, seq1 in enumerate(sequences):
            for j, seq2 in enumerate(sequences[i+1:], i+1):
                score, homology = self.pairwise_alignment(seq1, seq2)
                print(f"比对文件 {Path(file_paths[i]).name} 和 {Path(file_paths[j]).name}:\n比对得分: {score}\n同源性百分比: {homology:.2f}%\n")

    def process_command(self, args: argparse.Namespace) -> None:
        """处理序列相关命令

        Args:
            args: 命令行参数
        """
        if args.command == "fetch":
            try:
                with open(args.input_file, 'r', encoding='utf-8') as f:
                    genes = [line.strip() for line in f if line.strip()]
            except IOError as e:
                self.logger.log_error(f"读取输入文件失败: {args.input_file} - {str(e)}", self.config.error_log_path)
                return

            if genes:
                gene_uniprot_pairs = list(self.fetch_protein_sequences(genes, args.output_dir))
                if gene_uniprot_pairs:
                    self.fetch_domain_information(gene_uniprot_pairs, args.output_dir)
                else:
                    self.logger.log_error("未成功获取任何基因的UniProt ID，无法获取结构域信息。", self.config.error_log_path)
            else:
                self.logger.log_error("输入文件中未找到有效的基因名称。", self.config.error_log_path)

        elif args.command == "extract":
            if not os.path.isfile(args.fasta_file):
                self.logger.log_error(f"错误：文件 {args.fasta_file} 不存在。", self.config.error_log_path)
                return
            
            extracted_seq = self.extract_subsequence(args.fasta_file, args.start, args.end)
            if extracted_seq is None:
                return # Error already logged by extract_subsequence

            # Re-read header after successful extraction to ensure it's available
            header, _ = CommonUtils.read_fasta(args.fasta_file, return_header=True)
            if not isinstance(header, str):
                self.logger.log_error(f"无法从文件 {args.fasta_file} 获取头部信息。", self.config.error_log_path)
                return

            output_file = CommonUtils.generate_output_filename(args.fasta_file, ".fasta", args.start, args.end)
            CommonUtils.save_sequence(header, extracted_seq, output_file)
            print(f"已保存截取序列至 {output_file}")

        elif args.command == "mutate":
            if not os.path.isfile(args.fasta_file):
                self.logger.log_error(f"错误：文件 {args.fasta_file} 不存在。", self.config.error_log_path)
                return

            self.validate_input(args.pos, args.aa)
            try:
                record = next(SeqIO.parse(args.fasta_file, "fasta"))
            except Exception as e:
                self.logger.log_error(f"解析FASTA文件 {args.fasta_file} 失败: {str(e)}", self.config.error_log_path)
                return

            mutated_record = self.perform_mutations(record, args.pos, args.aa)
            mutation_info = "".join([f"{p}{a}" for p, a in zip(args.pos, args.aa)]) # Changed to p and a for clarity
            output_file = CommonUtils.generate_output_filename(os.path.splitext(args.fasta_file)[0], ".fasta", mutation_info)
            try:
                SeqIO.write(mutated_record, output_file, "fasta")
                print(f"已保存突变序列至 {output_file}")
            except IOError as e:
                self.logger.log_error(f"保存突变序列到文件失败: {output_file} - {str(e)}", self.config.error_log_path)

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
        try:
            config.read('exclude_residues.ini')
        except Exception as e:
            self.logger.log_error(f"读取 exclude_residues.ini 配置文件失败: {str(e)}", self.config.error_log_path)
            return []

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
            # Use ExPASy.get_sprot_raw for older SwissProt format, or UniProt API for modern data
            # Current implementation uses ExPASy, which relies on an older format.
            # A more robust approach might be to use the UniProt REST API for cross-references.
            # For now, keeping ExPASy but logging errors more consistently.
            handle = ExPASy.get_sprot_raw(uniprot_id)
            if handle is None:
                self.logger.log_error(f"无法从ExPASy获取UniProt ID {uniprot_id} 的原始数据。", self.config.error_log_path)
                return []
            record = SwissProt.read(handle)
            pdb_ids = []
            for cross_ref in record.cross_references:
                if cross_ref[0] == 'PDB':
                    pdb_ids.append(cross_ref[1])
            return pdb_ids
        except Exception as e:
            self.logger.log_error(f"获取UniProt ID {uniprot_id} 的PDB ID时出错: {e}", self.config.error_log_path)
            return []

    def download_pdb_files(self, pdb_ids: List[str], output_dir: str) -> None:
        """下载PDB文件

        Args:
            pdb_ids: PDB ID列表
            output_dir: 输出目录
        """
        output_dir_path = Path(output_dir)
        output_dir_path.mkdir(parents=True, exist_ok=True)
        pdbl = PDBList()
        failed_files: List[str] = [] # Track failed downloads

        def download_single_file(pdb_id: str) -> None:
            try:
                # PDBList.retrieve_pdb_file returns the local path if successful
                local_path = pdbl.retrieve_pdb_file(pdb_id, pdir=output_dir_path, file_format='mmCif')
                if local_path:
                    print(f"成功下载CIF文件: {pdb_id}.cif 到 {local_path}")
                else:
                    self.logger.log_error(f"下载CIF文件 {pdb_id} 失败，未返回路径。", self.config.error_log_path)
                    failed_files.append(pdb_id)
            except Exception as e:
                error_msg = f"下载CIF文件 {pdb_id} 时出错: {e}"
                self.logger.log_error(error_msg, self.config.error_log_path)
                failed_files.append(pdb_id)

        CommonUtils.parallel_executor(download_single_file, pdb_ids)

        if failed_files:
            failed_log_path = output_dir_path / 'download_failed.txt'
            with open(failed_log_path, 'a', encoding='utf-8') as f:
                for failed_file_id in failed_files:
                    f.write(f"{failed_file_id}.cif\n")
            print(f"部分CIF文件下载失败，详情请查看 {failed_log_path}")

    def extract_ligands_from_pdb(self, folder_path: str) -> Dict[str, List[str]]:
        """从.cif文件中提取配体信息

        Args:
            folder_path: 文件夹路径

        Returns:
            配体到PDB ID的映射字典
        """
        parser = MMCIFParser()
        ligand_pdb_dict: Dict[str, List[str]] = {}

        # Ensure folder_path exists and is a directory
        if not Path(folder_path).is_dir():
            self.logger.log_error(f"指定的PDB文件夹路径不存在或不是目录: {folder_path}", self.config.error_log_path)
            return {}

        for filename in os.listdir(folder_path):
            if filename.endswith('.cif'):
                pdb_id = filename.split('.')[0]
                file_path = os.path.join(folder_path, filename)
                try:
                    structure = parser.get_structure(pdb_id, file_path)
                    ligands: List[str] = []
                    for model in structure:
                        for chain in model:
                            for residue in chain:
                                # Check if it's a heteroatom residue (often starts with H_ in Biopython)
                                if residue.id[0].startswith('H_'):
                                    residue_name = residue.resname.strip()
                                    if residue_name not in self.exclude_residues and residue_name not in ligands:
                                        ligands.append(residue_name)
                except Exception as e:
                    error_msg = f"解析 {filename} 时出错: {e}"
                    self.logger.log_error(error_msg, self.config.error_log_path)
                    continue

                for ligand in ligands:
                    ligand_pdb_dict.setdefault(ligand, []).append(pdb_id)
        return ligand_pdb_dict

    @staticmethod
    def write_ligand_info_to_md(ligand_pdb_dict: Dict[str, List[str]], output_dir: str) -> None:
        """将配体信息写入pdb_ligand.md文件

        Args:
            ligand_pdb_dict: 配体到PDB ID的映射字典
            output_dir: 输出目录
        """
        output_file_path = Path(output_dir) / 'pdb_ligand.md'
        try:
            with open(output_file_path, 'w', encoding='utf-8') as outfile:
                outfile.write('| Ligands | PDB ID |\n')
                outfile.write('| --- | --- |\n')
                for ligand, pdb_ids in ligand_pdb_dict.items():
                    pdb_id_str = ', '.join(sorted(list(set(pdb_ids)))) # Use set to remove duplicates, then sort
                    outfile.write(f'| {ligand} | {pdb_id_str} |\n')
            print(f'结果已写入 {output_file_path}')
        except IOError as e:
            Logger.log_error(f"写入配体信息到MD文件失败: {output_file_path} - {str(e)}", ConfigManager().error_log_path)
            raise

    def move_files_based_on_ligands(self, folder_path: str, ligand_pdb_dict: Dict[str, List[str]]) -> None:
        """根据配体信息移动CIF文件到'no_ligand'或'with_ligands'子文件夹

        Args:
            folder_path: 包含CIF文件的文件夹路径
            ligand_pdb_dict: 配体到PDB ID的映射字典
        """
        no_ligand_dir = Path(folder_path) / 'no_ligand'
        with_ligands_dir = Path(folder_path) / 'with_ligands'

        no_ligand_dir.mkdir(exist_ok=True)
        with_ligands_dir.mkdir(exist_ok=True)

        pdb_ids_with_ligands = set()
        for pdb_ids in ligand_pdb_dict.values():
            pdb_ids_with_ligands.update(pdb_ids)

        for filename in os.listdir(folder_path):
            if filename.endswith('.cif'):
                pdb_id = filename.split('.')[0]
                src_path = Path(folder_path) / filename
                
                if pdb_id in pdb_ids_with_ligands:
                    dst_path = with_ligands_dir / filename
                    target_dir_name = 'with_ligands'
                else:
                    dst_path = no_ligand_dir / filename
                    target_dir_name = 'no_ligand'
                
                try:
                    # Only move if the file is not already in the target subfolder
                    if src_path.parent != dst_path.parent:
                        shutil.move(src_path, dst_path)
                        print(f"已将 {filename} 移动到 {target_dir_name}")
                except Exception as e:
                    self.logger.log_error(f"移动文件 {filename} 到 {target_dir_name} 时出错: {e}", self.config.error_log_path)

    def download_ligand_json(self, unique_ligands: Set[str], output_dir: str) -> None:
        """通过API并行查询配体的JSON文件并下载保存到json子文件夹

        Args:
            unique_ligands: 唯一配体集合
            output_dir: 输出目录
        """
        json_dir = Path(output_dir) / 'json'
        json_dir.mkdir(exist_ok=True)

        def download_single_ligand(ligand: str) -> None:
            retry_strategy = Retry(total=3, backoff_factor=1, status_forcelist=[500, 502, 503, 504])
            adapter = HTTPAdapter(max_retries=retry_strategy)
            session = requests.Session()
            session.mount('http://', adapter)
            session.mount('https://', adapter)

            try:
                url = f"https://data.rcsb.org/rest/v1/core/chemcomp/{ligand}"
                response = session.get(url, timeout=10)
                if response.status_code == 200:
                    data = response.json()
                    json_path = json_dir / f'{ligand}.json'
                    with open(json_path, 'w', encoding='utf-8') as json_file:
                        json.dump(data, json_file, ensure_ascii=False, indent=4)
                    print(f"JSON file for {ligand} has been saved to {json_path}")
                elif response.status_code == 404:
                    self.logger.log_error(f"配体 {ligand} 未在RCSB PDB找到 (404 Not Found)。", self.config.error_log_path)
                else:
                    self.logger.log_error(f"请求配体 {ligand} 失败，状态码: {response.status_code}", self.config.error_log_path)
            except requests.exceptions.RequestException as e:
                self.logger.log_error(f"下载配体 {ligand} 时发生网络错误: {e}", self.config.error_log_path)
            except Exception as e:
                self.logger.log_error(f"下载配体 {ligand} 时出错: {e}", self.config.error_log_path)
            finally:
                session.close()

        CommonUtils.parallel_executor(download_single_ligand, list(unique_ligands))

    def write_chemical_info_to_md(self, unique_ligands: Set[str], output_dir: str) -> None:
        """从json子文件夹中读取JSON文件并写入chemical_components_info.md

        Args:
            unique_ligands: 唯一配体集合
            output_dir: 输出目录
        """
        json_dir = Path(output_dir) / 'json'
        output_md_file = Path(output_dir) / 'chemical_components_info.md'
        
        try:
            with open(output_md_file, 'w', encoding='utf-8') as md_file:
                md_file.write("| Chemical Component ID | name | formula | formula_weight | Canonical Smiles |\n")
                md_file.write("| --- | --- | --- | --- | --- |\n")
                for ligand in unique_ligands:
                    json_path = json_dir / f'{ligand}.json'
                    try:
                        data = CommonUtils.load_from_json(str(json_path))
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
                        self.logger.log_error(f"JSON file for {ligand} not found in {json_dir}", self.config.error_log_path)
                    except Exception as e:
                        self.logger.log_error(f"读取或解析配体 {ligand} 的JSON文件失败: {str(e)}", self.config.error_log_path)
            print(f"MD file {output_md_file} has been created successfully.")
        except IOError as e:
            self.logger.log_error(f"创建或写入 chemical_components_info.md 文件失败: {str(e)}", self.config.error_log_path)
            raise

    def extract_ligands_coordinates(self, output_dir: str) -> None:
        """直接通过文本解析.cif文件提取配体坐标并保存成cif文件（并行版）

        Args:
            output_dir: 输出目录
        """
        with_ligands_dir = Path(output_dir) / 'with_ligands'
        ligands_dir = with_ligands_dir / 'ligands'
        ligands_dir.mkdir(exist_ok=True)

        cif_files = glob(str(with_ligands_dir / '*.cif'))
        if not cif_files:
            self.logger.log_error(f"在 {with_ligands_dir} 中未找到任何.cif文件进行配体坐标提取。", self.config.error_log_path)
            return

        def parse_cif_atom_sites(cif_content: str) -> Tuple[Dict[str, List[List[Any]]], Optional[str]]:
            atom_site_fields: List[str] = []
            atom_site_data: List[List[str]] = []
            in_loop = False
            data_block_id: Optional[str] = None

            lines = cif_content.split('\n')
            for line in lines:
                line = line.strip()
                if line.startswith('data_'):
                    data_block_id = line[5:] # Extract data block ID
                elif line.startswith('loop_'):
                    in_loop = True
                    atom_site_fields = [] # Reset fields for new loop
                    atom_site_data = [] # Reset data for new loop
                elif in_loop:
                    if line.startswith('_atom_site.'):
                        atom_site_fields.append(line.split('.')[1])
                    elif line and not line.startswith('#'):
                        # This assumes data lines immediately follow field definitions in the loop
                        # And that data lines are space-separated
                        atom_site_data.append(line.split())
                    elif not line and not atom_site_fields: # Empty line after loop_ or before first _atom_site.
                        pass # Continue looking for fields
                    elif line.startswith('#') or not line: # End of loop or end of data block
                        in_loop = False
                        if atom_site_fields and atom_site_data: # If we have collected fields and data for a loop
                            break # Assume we found the main atom_site loop, exit
                # If not in loop and not processing atom_site fields, just continue

            if not atom_site_fields or not atom_site_data:
                return {}, "未找到有效的_atom_site循环或数据"

            field_index = {field: idx for idx, field in enumerate(atom_site_fields)}
            required_fields = ['group_PDB', 'auth_comp_id', 'auth_asym_id', 'auth_seq_id', 
                              'Cartn_x', 'Cartn_y', 'Cartn_z', 'id', 'type_symbol']

            if not all(f in field_index for f in required_fields):
                return {}, f"缺少必要的_atom_site字段: {', '.join([f for f in required_fields if f not in field_index])}"

            ligands: Dict[str, List[List[Any]]] = {}
            for data_row in atom_site_data:
                try:
                    group_pdb = data_row[field_index['group_PDB']]
                    if group_pdb != 'HETATM':
                        continue

                    residue_name = data_row[field_index['auth_comp_id']].strip()
                    if residue_name in self.exclude_residues:
                        continue

                    x = float(data_row[field_index['Cartn_x']])
                    y = float(data_row[field_index['Cartn_y']])
                    z = float(data_row[field_index['Cartn_z']])

                    atom_fields = [
                        group_pdb,
                        residue_name,
                        data_row[field_index['auth_asym_id']],
                        data_row[field_index['auth_seq_id']],
                        x, y, z,
                        data_row[field_index['id']],
                        data_row[field_index['type_symbol']]
                    ]

                    ligands.setdefault(residue_name, []).append(atom_fields)
                except (ValueError, IndexError) as e:
                    self.logger.log_error(f"解析CIF行数据失败: {data_row} - {str(e)}", self.config.error_log_path)
                    continue
            return ligands, None

        def process_single_cif(cif_file: str) -> None:
            pdb_id = Path(cif_file).stem # Get filename without extension
            try:
                with open(cif_file, 'r', encoding='utf-8') as f:
                    cif_content = f.read()

                ligands, err = parse_cif_atom_sites(cif_content)
                if err:
                    raise ValueError(err)

                if not ligands:
                    # This case should ideally not happen if file was moved to 'with_ligands'
                    # unless it contains only excluded residues or no HETATM records.
                    self.logger.log_error(f"文件 {pdb_id}.cif 中未找到可提取的配体原子坐标。", self.config.error_log_path)
                    return

                for residue_name, atom_fields_list in ligands.items():
                    ligand_cif_file = ligands_dir / f'{pdb_id}_{residue_name}.cif'

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
                        # Format coordinates to 3 decimal places
                        cif_lines.append(
                            f"{atom_fields[0]} {atom_fields[1]} {atom_fields[2]} {atom_fields[3]} "
                            f"{atom_fields[4]:.3f} {atom_fields[5]:.3f} {atom_fields[6]:.3f} "
                            f"{atom_fields[7]} {atom_fields[8]}\n"
                        )
                    cif_lines.append("#\n") # End of data block

                    with open(ligand_cif_file, 'w', encoding='utf-8') as f:
                        f.writelines(cif_lines)
                    print(f"成功保存配体 {residue_name} 的CIF文件: {ligand_cif_file}")

            except Exception as e:
                error_msg = f"处理CIF文件 {cif_file} 时出错: {e}"
                self.logger.log_error(error_msg, self.config.error_log_path)

        CommonUtils.parallel_executor(process_single_cif, cif_files)

    def process(self) -> None:
        """处理PDB相关任务"""
        for section in self.config.pdb_config.sections():
            output_dir = self.config.pdb_config.get(section, 'output', fallback=None)
            uniprot_id = self.config.pdb_config.get(section, 'uniprot', fallback=None)

            if not output_dir:
                self.logger.log_error(f"PDB配置文件中 {section} 部分缺少 'output' 路径。", self.config.error_log_path)
                continue
            if not uniprot_id:
                self.logger.log_error(f"PDB配置文件中 {section} 部分缺少 'uniprot' ID。", self.config.error_log_path)
                continue
            
            print(f"\n--- 开始处理 PDB 任务: {section} (UniProt ID: {uniprot_id}) ---")
            
            pdb_ids = self.get_pdb_ids_from_uniprot(uniprot_id)
            if pdb_ids:
                print(f"为 UniProt ID {uniprot_id} 找到 PDB IDs: {', '.join(pdb_ids)}")
                self.download_pdb_files(pdb_ids, output_dir)
            else:
                self.logger.log_error(f"任务 {section}: 未找到对应的PDB ID。", self.config.error_log_path)
                continue # Skip further processing for this section if no PDB IDs found

            # Re-check if any CIF files were downloaded before proceeding
            downloaded_cif_files = glob(str(Path(output_dir) / '*.cif'))
            if not downloaded_cif_files:
                self.logger.log_error(f"在 {output_dir} 中未找到下载的CIF文件，跳过后续配体处理。", self.config.error_log_path)
                continue

            ligand_pdb_dict = self.extract_ligands_from_pdb(output_dir)
            
            # Moved file organization before writing ligand info to MD, as files are moved *out* of output_dir
            # This ensures only relevant files are in 'output_dir' for ligand extraction later.
            # However, the previous logic was to move *after* extraction, which means extract_ligands_from_pdb
            # should operate on the initial download directory.
            # Let's adjust: extract, then write MD, then move. This maintains the flow.
            
            if ligand_pdb_dict:
                self.write_ligand_info_to_md(ligand_pdb_dict, output_dir)
            else:
                print(f"在 {output_dir} 中未找到任何配体信息。")
            
            # Move files after all extractions from the original directory are done
            self.move_files_based_on_ligands(output_dir, ligand_pdb_dict)
            
            # Read ligands from the generated MD file to ensure consistency
            md_file_path = Path(output_dir) / 'pdb_ligand.md'
            unique_ligands: Set[str] = set()
            if md_file_path.exists():
                try:
                    with open(md_file_path, 'r', encoding='utf-8') as f:
                        lines = f.readlines()
                    for line in lines[2:]: # Skip header and separator
                        parts = [p.strip() for p in line.strip().split('|') if p.strip()]
                        if len(parts) >= 1:
                            unique_ligands.add(parts[0])
                except IOError as e:
                    self.logger.log_error(f"读取 {md_file_path} 失败: {str(e)}", self.config.error_log_path)
            else:
                self.logger.log_error(f"配体信息文件 {md_file_path} 不存在，无法获取唯一配体列表。", self.config.error_log_path)

            if unique_ligands:
                print(f"找到唯一配体: {', '.join(sorted(list(unique_ligands)))}")
                self.download_ligand_json(unique_ligands, output_dir)
                self.write_chemical_info_to_md(unique_ligands, output_dir)
                self.extract_ligands_coordinates(output_dir)
            else:
                print("未找到任何配体，跳过配体JSON下载和坐标提取。")

        print('\n所有PDB任务完成。')


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
        self.session = self._setup_session()

    def _setup_session(self) -> requests.Session:
        """设置带有重试策略的requests会话"""
        retry_strategy = Retry(
            total=5,
            backoff_factor=1,
            status_forcelist=[429, 500, 502, 503, 504],
            allowed_methods=["HEAD", "GET", "OPTIONS"]
        )
        adapter = HTTPAdapter(max_retries=retry_strategy)
        session = requests.Session()
        session.mount('http://', adapter)
        session.mount('https://', adapter)
        session.headers.update({
            'Accept': 'application/json',
            'User-Agent': 'PyPDA/1.0 (https://github.com/example/pypda; pypda@example.com)' # Good practice to provide contact info
        })
        return session

    def get_uniprot_data(self, accession: str) -> Optional[Dict[str, Any]]:
        """根据UniProt编号从API获取数据

        Args:
            accession: UniProt蛋白质编号

        Returns:
            如果成功获取，返回蛋白质信息的字典；否则返回None。
        """
        url = f'{self.api_base_url}/{accession}'
        try:
            print(f"正在从 {url} 获取数据...")
            response = self.session.get(url, timeout=30) # Increased timeout
            response.raise_for_status() # Raise HTTPError for bad responses (4xx or 5xx)
            return response.json()
        except requests.exceptions.HTTPError as e:
            self.logger.log_error(f"HTTP请求错误 (UniProt API - {url}): {e.response.status_code} - {e.response.text}", self.config.error_log_path)
        except requests.exceptions.ConnectionError:
            self.logger.log_error(f"网络连接错误 (UniProt API - {url})", self.config.error_log_path)
        except requests.exceptions.Timeout:
            self.logger.log_error(f"请求超时 (UniProt API - {url})", self.config.error_log_path)
        except json.JSONDecodeError:
            self.logger.log_error(f"UniProt API返回无效JSON (UniProt API - {url}): {response.text[:200]}...", self.config.error_log_path)
        except Exception as e:
            self.logger.log_error(f"获取UniProt数据失败: {str(e)} (UniProt API - {url})", self.config.error_log_path)
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
                         'TISSUE SPECIFICITY', 'SUBCELLULAR LOCATION', 'PTM', 'SIMILARITY', 'DISEASE', 'INTERACTION'] # Added DISEASE, INTERACTION for completeness
        # Create a case-insensitive mapping dictionary
        comment_type_map = {ct.lower(): ct for ct in comment_types}
        info['comments'] = {ct: [] for ct in comment_types}
        
        for comment in comments:
            # Convert API's commentType to lowercase for matching
            ct_lower = comment.get('commentType', '').lower()
            if ct_lower in comment_type_map:
                ct = comment_type_map[ct_lower]
                info['comments'][ct].append(comment)
            elif 'INTERACTION' in comment_type_map and ct_lower == 'interaction': # Handle explicit interaction
                 info['comments']['INTERACTION'].append(comment)
            # else:
            #     print(f"Unknown comment type: {comment.get('commentType')}") # For debugging new types

        # 5. 蛋白质特征
        features = data.get('features', [])
        info['features'] = {
            'counts': {},
            'detailed': {}
        }
        detailed_feature_types = ["Chain", "Region", "Active site", "Binding site", "Modified residue", "Mutagenesis", "Domain"] # Added Domain
        for dt in detailed_feature_types:
            info['features']['detailed'][dt] = []

        for feature in features:
            ft = feature.get('type')
            info['features']['counts'][ft] = info['features']['counts'].get(ft, 0) + 1
            if ft in detailed_feature_types:
                info['features']['detailed'][ft].append(feature)

        # 6. 蛋白质相互作用 (separate from comments for clarity if needed, but comments already cover it)
        # The 'interactions' key in the top-level data provides a more structured view
        info['interactions'] = data.get('interactions', [])

        # 7. 关键词
        info['keywords'] = [kw.get('name') for kw in data.get('keywords', []) if kw.get('name')]

        # 8. 参考文献
        info['references'] = data.get('references', [])

        # 9. 交叉引用
        cross_references = data.get('uniProtKBCrossReferences', [])
        important_databases = ["PDB", "DrugBank", "GO", "Reactome", "HGNC", "GeneID", "KEGG", "AlphaFoldDB", "IntAct"] # Added AlphaFoldDB, IntAct
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
                f.write(f"# UniProt 蛋白质信息分析报告 - {info['basic_info'].get('primaryAccession', 'N/A')}\n\n")

                # 1. 基本识别信息
                f.write("## 1. 基本识别信息\n")
                bi = info['basic_info']
                f.write(f"- **条目类型**: {bi.get('entryType', 'N/A')}\n")
                f.write(f"- **主要登录号**: {bi.get('primaryAccession', 'N/A')}\n")
                if bi.get('secondaryAccessions'):
                    f.write(f"- **次要登录号**: {', '.join(bi['secondaryAccessions'])}\n")
                f.write(f"- **UniProtKB ID**: {bi.get('uniProtkbId', 'N/A')}\n")
                f.write(f"- **注释评分**: {bi.get('annotationScore', 'N/A')}\n")
                audit = bi.get('entryAudit', {})
                f.write(f"- **首次公开日期**: {audit.get('firstPublicDate', 'N/A')}\n")
                f.write(f"- **最后注释更新日期**: {audit.get('lastAnnotationUpdateDate', 'N/A')}\n")
                f.write(f"- **最后序列更新日期**: {audit.get('lastSequenceUpdateDate', 'N/A')}\n")
                f.write(f"- **条目版本**: {audit.get('entryVersion', 'N/A')}\n")
                f.write(f"- **序列版本**: {audit.get('sequenceVersion', 'N/A')}\n\n")

                # 2. 生物学背景信息
                f.write("## 2. 生物学背景信息\n")
                bio = info['biology_info']
                f.write(f"- **科学名称**: {bio.get('scientificName', 'N/A')}\n")
                f.write(f"- **常用名称**: {bio.get('commonName', 'N/A')}\n")
                f.write(f"- **分类ID**: {bio.get('taxonId', 'N/A')}\n")
                if bio.get('lineage'):
                    f.write(f"- **生物学谱系**: {' -> '.join(bio['lineage'])}\n")
                f.write(f"- **蛋白质存在证据**: {bio.get('proteinExistence', 'N/A')}\n\n")

                # 3. 蛋白质描述和基因信息
                f.write("## 3. 蛋白质描述和基因信息\n")
                pd = info['protein_desc']
                rn = pd['recommendedName']
                f.write(f"- **推荐全名**: {rn.get('fullName', 'N/A')}\n")
                if rn.get('shortNames'):
                    f.write(f"- **推荐简称**: {', '.join(rn['shortNames'])}\n")
                if rn.get('ecNumbers'):
                    f.write(f"- **EC 编号**: {', '.join(rn['ecNumbers'])}\n")

                alternative_names = pd.get('alternativeNames', [])
                if alternative_names:
                    f.write("- **备选名称**:\n")
                    for alt_name in alternative_names:
                        alt_full_name = alt_name.get('fullName', {}).get('value', 'N/A')
                        alt_ec_numbers = [ec.get('value') for ec in alt_name.get('ecNumbers', []) if ec.get('value')]
                        f.write(f"  - {alt_full_name}" + (f" (EC: {', '.join(alt_ec_numbers)})" if alt_ec_numbers else "") + "\n")

                genes = pd.get('genes', [])
                if genes:
                    for gene in genes:
                        gene_name = gene.get('geneName', {}).get('value', 'N/A')
                        f.write(f"- **基因名称**: {gene_name}\n")
                        gene_synonyms = [syn.get('value') for syn in gene.get('synonyms', []) if syn.get('value')]
                        if gene_synonyms:
                            f.write(f"  - **基因同义词**: {', '.join(gene_synonyms)}\n")
                f.write("\n")

                # 4. 功能和活性注释
                f.write("## 4. 功能和活性注释\n")
                comments_data = info['comments']
                comment_order = ['FUNCTION', 'CATALYTIC ACTIVITY', 'COFACTOR', 'ACTIVITY REGULATION',
                                 'TISSUE SPECIFICITY', 'SUBCELLULAR LOCATION', 'PTM', 'SIMILARITY',
                                 'DISEASE', 'INTERACTION'] # Ensure consistent order
                
                found_any_comment = False
                for ct in comment_order:
                    if comments_data.get(ct):
                        found_any_comment = True
                        f.write(f"### {ct.replace('_', ' ').title()}\n") # Format to Title Case
                        for comment in comments_data[ct]:
                            texts = []
                            # Handle different structures of 'text' field (list of dicts or single dict/string)
                            if 'texts' in comment and isinstance(comment['texts'], list):
                                for text_item in comment['texts']:
                                    if isinstance(text_item, dict) and 'value' in text_item:
                                        texts.append(text_item['value'])
                                    elif isinstance(text_item, str):
                                        texts.append(text_item)
                            elif 'text' in comment: # For comments with a single 'text' field (e.g., PTM)
                                if isinstance(comment['text'], dict) and 'value' in comment['text']:
                                    texts.append(comment['text']['value'])
                                elif isinstance(comment['text'], str):
                                    texts.append(comment['text'])

                            if texts:
                                f.write(f"- {'; '.join(texts)}\n")
                            
                            # Additional details for specific comment types
                            if ct == 'CATALYTIC ACTIVITY':
                                reaction = comment.get('reaction')
                                if reaction:
                                    f.write(f"  - **Reaction**: {reaction.get('name', 'N/A')}\n")
                            elif ct == 'DISEASE':
                                disease_name = comment.get('disease', {}).get('diseaseName', 'N/A')
                                acronym = comment.get('disease', {}).get('acronym', 'N/A')
                                f.write(f"  - **Disease**: {disease_name} ({acronym})\n")
                            elif ct == 'SUBCELLULAR LOCATION':
                                for location in comment.get('locations', []):
                                    f.write(f"  - **Location**: {location.get('location', {}).get('value', 'N/A')}\n")
                                    for topology in location.get('topologies', []):
                                        f.write(f"    - **Topology**: {topology.get('value', 'N/A')}\n")
                                    for orientation in location.get('orientations', []):
                                        f.write(f"    - **Orientation**: {orientation.get('value', 'N/A')}\n")

                if not found_any_comment:
                    f.write("- 未找到功能和活性注释信息。\n")
                f.write("\n")

                # 5. 蛋白质特征
                f.write("## 5. 蛋白质特征\n")
                features = info['features']
                f.write("### 特征类型统计\n")
                if features['counts']:
                    for ft, count in features['counts'].items():
                        f.write(f"- {ft}: {count}\n")
                else:
                    f.write("- 未找到特征类型统计信息。\n")

                f.write("\n### 详细特征信息\n")
                found_any_detailed_feature = False
                for dt in features['detailed']:
                    if features['detailed'][dt]:
                        found_any_detailed_feature = True
                        f.write(f"#### {dt}\n")
                        for feature in features['detailed'][dt]:
                            desc = feature.get('description', 'N/A')
                            if isinstance(desc, dict): # Sometimes description can be a dict with 'value'
                                desc = desc.get('value', 'N/A')

                            loc = feature.get('location', {})
                            begin = loc.get('start', {}).get('value', 'N/A')
                            end = loc.get('end', {}).get('value', 'N/A')
                            
                            feature_line = f"- {desc} (位置: {begin}-{end})"
                            # Add original amino acid for 'Modified residue' or 'Mutagenesis'
                            if dt == 'Modified residue' or dt == 'Mutagenesis':
                                original_aa = feature.get('original', 'N/A')
                                feature_line += f", 原始氨基酸: {original_aa}"
                            
                            f.write(feature_line + "\n")
                if not found_any_detailed_feature:
                    f.write("- 未找到详细特征信息。\n")
                f.write("\n")

                # 6. 蛋白质相互作用
                f.write("## 6. 蛋白质相互作用\n")
                interactions = info['interactions']
                if interactions:
                    for i, interaction in enumerate(interactions, 1):
                        interactor = interaction.get('interactor', {})
                        interactor_id = interactor.get('uniprotId', 'N/A')
                        interactor_name = interactor.get('geneDisplayName', interactor.get('name', 'N/A')) # Use geneDisplayName if available
                        experiments = interaction.get('experiments', 0)
                        methods = [m.get('name', 'N/A') for m in interaction.get('methods', [])]

                        f.write(f"- 相互作用蛋白 {i}: **{interactor_name}** (UniProt ID: {interactor_id})\n")
                        f.write(f"  - 实验数量: {experiments}\n")
                        if methods:
                            f.write(f"  - 相互作用方法: {', '.join(methods)}\n")
                else:
                    f.write("- 未找到相互作用信息\n")
                f.write("\n")

                # 7. 关键词
                f.write("## 7. 关键词\n")
                keywords = info['keywords']
                if keywords:
                    f.write(f"- {', '.join(sorted(keywords))}\n") # Sort keywords for consistent output
                else:
                    f.write("- 未找到关键词信息\n")
                f.write("\n")

                # 8. 参考文献
                f.write("## 8. 参考文献\n")
                references = info['references']
                if references:
                    for i, ref in enumerate(references[:10], 1): # Display up to 10 references
                        citation = ref.get('citation', {})
                        authors_list = citation.get('authors', [])
                        authors = ", ".join(authors_list) if authors_list else 'N/A'
                        title = citation.get('title', 'N/A')
                        journal = citation.get('journal', 'N/A')
                        year = citation.get('publicationDate', 'N/A')[:4]  # Extract year
                        pubmed_id = ""
                        for xref in citation.get('uniProtKBCrossReferences', []):
                            if xref.get('database') == 'PubMed':
                                pubmed_id = f" [PMID:{xref.get('id', '')}]"
                                break
                        f.write(f"- [{i}] {authors}, \"{title}\", *{journal}*, {year}{pubmed_id}\n")
                    if len(references) > 10:
                        f.write(f"- 显示前10篇，共{len(references)}篇参考文献\n")
                else:
                    f.write("- 未找到参考文献信息\n")
                f.write("\n")

                # 9. 交叉引用
                f.write("## 9. 交叉引用\n")
                cross_refs = info['cross_references']
                found_any_xref = False
                for db in sorted(cross_refs.keys()): # Sort databases for consistent output
                    if cross_refs[db]:
                        found_any_xref = True
                        f.write(f"### {db}\n")
                        ids = [xref.get('id') for xref in cross_refs[db] if xref.get('id')]
                        # Limit displayed IDs to avoid overly long lines
                        display_ids = ids[:10]
                        f.write(f"- {', '.join(display_ids)}")
                        if len(ids) > 10:
                            f.write(f" ... (共{len(ids)}个条目)")
                        f.write("\n")
                if not found_any_xref:
                    f.write("- 未找到交叉引用信息。\n")
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
                else:
                    f.write("- 未找到氨基酸序列信息。\n")

            print(f"Markdown报告已生成至 {md_filename}")
        except IOError as e:
            Logger.log_error(f"生成Markdown报告失败: {md_filename} - {str(e)}", ConfigManager().error_log_path)
            raise


class PypdaApp:
    """Pypda应用主类，负责命令行参数解析和工具调度"""
    def __init__(self):
        self.config = ConfigManager()
        self.logger = Logger()
        self.seq_processor = SequenceProcessor(self.config, self.logger)
        self.pdb_processor = PDBProcessor(self.config, self.logger)
        self.uniprot_api = UniProtAPI(self.config, self.logger)

    def _setup_seq_parser(self, subparsers: argparse._SubParsersAction) -> None:
        """设置序列分析工具的子命令解析器"""
        seq_parser = subparsers.add_parser(
            "seq", 
            help="蛋白质序列处理工具，支持序列获取、提取、突变和比对。",
            formatter_class=argparse.RawTextHelpFormatter # For better help formatting
        )
        seq_subparsers = seq_parser.add_subparsers(dest="command", required=True, help="序列处理命令")

        # seq fetch命令
        fetch_parser = seq_subparsers.add_parser(
            "fetch", 
            help="批量获取蛋白质序列和结构域信息。",
            description="""
            从UniProt获取指定基因的人类蛋白质全长序列，并提取结构域信息。
            输入文件应包含每行一个HGNC基因名称。
            """
        )
        fetch_parser.add_argument("input_file", type=str, help="包含HGNC基因名称的文本文件路径，每行一个基因名称。")
        fetch_parser.add_argument(
            "output_dir", 
            type=str,
            default="protein_sequences", 
            nargs='?', # Make default an option
            help="输出目录，用于保存FASTA序列文件和结构域信息报告 (默认: protein_sequences)。"
        )

        # seq extract命令
        extract_parser = seq_subparsers.add_parser(
            "extract", 
            help="从FASTA文件中提取指定位置的子序列。",
            description="""
            从给定的FASTA文件（通常是蛋白质序列）中，根据起始和结束位置提取子序列。
            起始和结束位置均为1-based索引。
            """
        )
        extract_parser.add_argument("fasta_file", type=str, help="输入FASTA文件路径。")
        extract_parser.add_argument("start", type=int, help="子序列的起始位置 (1-based)。")
        extract_parser.add_argument("end", type=int, help="子序列的结束位置 (1-based)。")

        # seq mutate命令
        mut_parser = seq_subparsers.add_parser(
            "mutate", 
            help="对蛋白质序列执行点突变。",
            description="""
            对FASTA文件中的蛋白质序列执行一个或多个点突变。
            突变位置和新氨基酸列表必须一一对应。
            例如：--pos 10 20 --aa A G 表示将第10位突变为丙氨酸，第20位突变为甘氨酸。
            """
        )
        mut_parser.add_argument("fasta_file", type=str, help="输入FASTA文件路径。")
        mut_parser.add_argument(
            "--pos", 
            nargs="+", 
            type=int, 
            required=True, 
            help="一个或多个突变位置 (1-based)，用空格分隔。"
        )
        mut_parser.add_argument(
            "--aa", 
            nargs="+", 
            type=str, 
            required=True, 
            help="与突变位置对应的新的氨基酸单字母代码，用空格分隔。"
        )

        # seq align命令
        align_parser = seq_subparsers.add_parser(
            "align", 
            help="对多个蛋白质序列进行两两比对。",
            description="""
            对提供的所有FASTA文件中的蛋白质序列进行两两全局比对，
            并报告比对得分和同源性百分比。
            """
        )
        align_parser.add_argument("fasta_files", nargs="+", type=str, help="一个或多个FASTA文件路径，用于比对。")

    def _setup_pdb_parser(self, subparsers: argparse._SubParsersAction) -> None:
        """设置PDB文件处理工具的子命令解析器"""
        pdb_parser = subparsers.add_parser(
            "pdb", 
            help="PDB文件处理工具，支持PDB文件下载、配体信息提取和文件整理。",
            description="""
            根据配置文件 (pdb_config.ini) 中的UniProt ID，自动下载相关PDB文件 (mmCIF格式)，
            提取蛋白质中的配体信息，生成报告，并根据是否含有配体将PDB文件分类。
            同时，下载配体化学信息并提取配体坐标。
            """
        )
        # PDB command doesn't have subcommands, it just triggers the PDBProcessor.process()
        # No additional arguments are directly exposed via CLI for this command,
        # as its behavior is driven by pdb_config.ini.
        # Adding a dummy argument or making it a direct command without sub-subparsers
        # is a design choice. Given the current `process` method, no direct args are needed.

    def _setup_uniprot_parser(self, subparsers: argparse._SubParsersAction) -> None:
        """设置UniProt数据处理工具的子命令解析器"""
        uniprot_parser = subparsers.add_parser(
            "uniprot", 
            help="UniProt数据处理工具，支持从UniProt API获取数据及分析本地JSON文件。",
            formatter_class=argparse.RawTextHelpFormatter
        )
        uniprot_subparsers = uniprot_parser.add_subparsers(dest="command", required=True, help="UniProt数据命令")

        # uniprot fetch命令
        uniprot_fetch_parser = uniprot_subparsers.add_parser(
            "fetch", 
            help="从UniProt API获取蛋白质数据并生成详细报告。",
            description="""
            根据UniProt蛋白质编号，从UniProt REST API获取完整的蛋白质信息，
            并将其保存为JSON文件和Markdown格式的分析报告。
            """
        )
        uniprot_fetch_parser.add_argument("accession", type=str, help="UniProt蛋白质编号，例如: Q13547。")
        uniprot_fetch_parser.add_argument(
            "-o", "--output_dir", 
            type=str, 
            default="uniprot_reports", 
            help="保存JSON和Markdown报告的输出目录 (默认: uniprot_reports)。"
        )

        # uniprot analyze命令
        uniprot_analyze_parser = uniprot_subparsers.add_parser(
            "analyze", 
            help="分析现有UniProt JSON文件并生成报告。",
            description="""
            加载本地已有的UniProt蛋白质信息JSON文件，
            提取关键数据并生成Markdown格式的分析报告。
            """
        )
        uniprot_analyze_parser.add_argument(
            "-f", "--file", 
            type=str, 
            required=True, 
            help="要分析的UniProt蛋白质信息JSON文件路径。"
        )

    def setup_parser(self) -> argparse.ArgumentParser:
        """设置命令行参数解析器

        Returns:
            配置好的参数解析器
        """
        parser = argparse.ArgumentParser(
            description="蛋白质数据分析综合工具 (PyPDA)",
            formatter_class=argparse.RawTextHelpFormatter # For general help formatting
        )
        
        # Add a version argument
        parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0.0')

        subparsers = parser.add_subparsers(
            dest="tool", 
            required=True, 
            help="选择要使用的工具：序列处理 (seq), PDB文件处理 (pdb), 或UniProt数据分析 (uniprot)。"
        )

        self._setup_seq_parser(subparsers)
        self._setup_pdb_parser(subparsers)
        self._setup_uniprot_parser(subparsers)

        return parser

    def run(self) -> None:
        """运行应用程序"""
        parser = self.setup_parser()
        args = parser.parse_args()

        try:
            if args.tool == "seq":
                self.seq_processor.process_command(args)
            elif args.tool == "pdb":
                self.pdb_processor.process()
            elif args.tool == "uniprot":
                if args.command == "fetch":
                    # Create output directory for uniprot reports
                    output_dir_path = Path(args.output_dir)
                    output_dir_path.mkdir(parents=True, exist_ok=True)

                    data = self.uniprot_api.get_uniprot_data(args.accession)
                    if not data:
                        self.logger.log_error("无法获取UniProt数据。", self.config.error_log_path)
                        return
                    json_filename = CommonUtils.save_to_json(data, args.accession, str(output_dir_path))
                    analyzer = ProteinAnalyzer()
                    protein_info = analyzer.extract_protein_info(data)
                    md_filename = str(output_dir_path / (Path(json_filename).stem + '.md'))
                    ReportGenerator.generate_md_report(protein_info, md_filename)
                elif args.command == "analyze":
                    data = CommonUtils.load_from_json(args.file)
                    analyzer = ProteinAnalyzer()
                    protein_info = analyzer.extract_protein_info(data)
                    md_filename = str(Path(args.file).with_suffix('.md')) # Save MD in same directory as JSON
                    ReportGenerator.generate_md_report(protein_info, md_filename)
        except Exception as e:
            self.logger.log_error(f"应用程序运行过程中发生未捕获的错误: {str(e)}", self.config.error_log_path)
            # Optionally, re-raise for debugging during development: raise

if __name__ == "__main__":
    app = PypdaApp()
    app.run()
