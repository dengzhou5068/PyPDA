# 标准库导入
import argparse
import configparser
import io
import json
import os
import requests
import shutil
import sys
import warnings
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
from functools import lru_cache
from pathlib import Path
from typing import List, Dict, Tuple, Generator, Optional, Set, Any, TypedDict, Union

# 第三方库导入
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from Bio.PDB.PDBList import PDBList
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from requests.adapters import HTTPAdapter
from tqdm import tqdm
from urllib3.util.retry import Retry

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, DataStructs
except ImportError:
    warnings.warn("RDKit库未正确安装，将无法使用分子结构相似性计算功能。")

# 过滤 PDB 构建警告
warnings.filterwarnings("ignore", category=PDBConstructionWarning)


class ConfigManager:
    """配置文件管理类，负责读取和管理配置信息"""
    def __init__(self):
        self.ns = {'uniprot': 'http://uniprot.org/uniprot'}
        self.error_log_path = os.path.join(os.getcwd(), 'error.txt')
        self.info_log_path = os.path.join(os.getcwd(), 'info.txt')
        self.uniprot_api_base_url = 'https://rest.uniprot.org/uniprotkb/'


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
        print(f"ERROR: {error_msg}", file=sys.stderr) # Print to stderr for errors
        
    @staticmethod
    def log_info(info_msg: str, info_file_path: str) -> None:
        """记录信息日志到日志文件并打印

        Args:
            info_msg: 信息内容
            info_file_path: 信息日志文件路径
        """
        with open(info_file_path, 'a', encoding='utf-8') as f:
            f.write(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - INFO: {info_msg}\n")
        print(f"INFO: {info_msg}") # Print to stdout for info messages


class CommonUtils:
    """通用工具类，提供项目中常用的工具函数"""
    @staticmethod
    def read_fasta(file_path: Union[str, Path], return_header: bool = False) -> Union[str, Tuple[str, str]]:
        """读取FASTA文件，可选择返回头部信息

        Args:
            file_path: 文件路径
            return_header: 是否返回头部信息

        Returns:
            序列字符串，或(头部, 序列)元组
        """
        try:
            record = SeqIO.read(str(file_path), "fasta")
            if return_header:
                return record.description, str(record.seq)
            return str(record.seq)
        except FileNotFoundError:
            Logger.log_error(f"FASTA文件未找到: {file_path}", ConfigManager().error_log_path)
            raise
        except Exception as e:
            Logger.log_error(f"读取FASTA文件失败: {file_path} - {str(e)}", ConfigManager().error_log_path)
            raise

    @staticmethod
    def generate_output_filename(base_name: Union[str, Path], ext: str, *args: Any) -> str:
        """统一生成输出文件名

        Args:
            base_name: 基础文件名
            ext: 文件扩展名（含.）
            *args: 可变参数，用于生成附加信息

        Returns:
            生成的文件名
        """
        base_name_path = Path(base_name)
        if not args:
            return f"{base_name_path.stem}{ext}"
        info = "_".join(map(str, args)) if len(args) > 1 else str(args[0])
        return f"{base_name_path.stem}_{info}{ext}"

    @staticmethod
    def save_sequence(header: str, sequence: str, output_file: Union[str, Path]) -> None:
        """保存序列到FASTA文件

        Args:
            header: 序列头部信息
            sequence: 序列内容
            output_file: 输出文件路径
        """
        try:
            with open(output_file, 'w', encoding="utf-8") as file:
                # 确保头部信息以>符号开头
                if not header.startswith('>'):
                    header = '>' + header
                file.write(header + '\n')
                for i in range(0, len(sequence), 60):
                    file.write(sequence[i:i+60] + '\n')
            print(f"序列已保存至 {output_file}")
        except IOError as e:
            Logger.log_error(f"保存序列到文件失败: {output_file} - {str(e)}", ConfigManager().error_log_path)
            raise

    @staticmethod
    def parallel_executor(func: callable, items: List[Any], max_workers: Optional[int] = None, description: str = "Processing") -> None:
        """并行执行函数，并显示进度条

        Args:
            func: 要并行执行的函数
            items: 迭代参数列表
            max_workers: 最大工作线程数
            description: 进度条描述
        """
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            list(tqdm(executor.map(func, items), total=len(items), desc=description))

    @staticmethod
    def save_to_json(data: Dict[str, Any], accession: str, output_dir: Union[str, Path] = '.') -> Path:
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
            return json_filename
        except IOError as e:
            Logger.log_error(f"保存JSON文件失败: {json_filename} - {str(e)}", ConfigManager().error_log_path)
            raise

    @staticmethod
    def load_from_json(file_path: Union[str, Path]) -> Dict[str, Any]:
        """从JSON文件加载数据

        Args:
            file_path: JSON文件路径

        Returns:
            从文件加载的字典数据
        """
        if not Path(file_path).exists():
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

    @staticmethod
    def get_output_dir(base_dir: str, tool_name: str, input_name: str = "output", timestamp: str = None) -> str:
        """根据工具名称和输入名称生成带时间戳的输出目录

        Args:
            base_dir: 基础目录，如 "result"
            tool_name: 工具名称，如 "protein_sequences", "pdb_output", "uniprot_reports"
            input_name: 输入名称，如基因名或蛋白质名
            timestamp: 时间戳，如果为None则使用当前时间

        Returns:
            完整的输出目录路径
        """
        if timestamp is None:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_dir = Path(base_dir) / tool_name / f"{input_name}_{timestamp}"
        output_dir.mkdir(parents=True, exist_ok=True)
        return str(output_dir), timestamp


class SequenceProcessor:
    """序列处理类，负责蛋白质序列相关操作"""
    def __init__(self, config: ConfigManager, logger: Logger, uniprot_api: 'UniProtAPI'):
        self.config = config
        self.logger = logger
        self.uniprot_api = uniprot_api
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

        for gene in tqdm(genes, desc="Fetching sequences"):
            params = {
                "query": f"gene_exact:{gene} AND organism_id:9606",
                "format": "fasta",
                "fields": "accession,sequence"
            }
            try:
                response = self.uniprot_api.session.get(
                    f"{self.config.uniprot_api_base_url}search",
                    params=params,
                    timeout=30
                )
                response.raise_for_status()

                if not response.text.strip():
                    self.logger.log_error(f"未找到基因 {gene} 对应的蛋白质序列或响应为空", self.config.error_log_path)
                    continue

                output_file = output_dir_path / f"{gene}.fasta"
                try:
                    first_seq_record = next(SeqIO.parse(io.StringIO(response.text.strip()), "fasta"))
                    SeqIO.write(first_seq_record, output_file, "fasta")
                except StopIteration:
                    self.logger.log_error(f"UniProt API返回无效FASTA格式，针对基因 {gene}", self.config.error_log_path)
                    continue

                print(f"已将 {gene} 的第一个全长蛋白质序列保存到 {output_file}")
                uniprot_id = first_seq_record.id.split('|')[1] if '|' in first_seq_record.id else first_seq_record.id
                yield gene, uniprot_id

            except requests.exceptions.RequestException as e:
                self.logger.log_error(f"获取 {gene} 蛋白质序列失败: {e}", self.config.error_log_path)
            except Exception as e:
                self.logger.log_error(f"处理基因 {gene} 时发生未知错误: {e}", self.config.error_log_path)


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

        for gene, uniprot_id in tqdm(gene_uniprot_pairs, desc="Fetching domain info"):
            try:
                fasta_file = output_dir_path / f"{gene}.fasta"
                total_amino_acids: Union[int, str]
                if fasta_file.exists():
                    record = SeqIO.read(fasta_file, "fasta")
                    total_amino_acids = len(record.seq)
                else:
                    total_amino_acids = "未知"

                data = self.uniprot_api.get_uniprot_data(uniprot_id)
                if not data:
                    self.logger.log_error(f"无法获取UniProt ID {uniprot_id} 的详细数据，跳过结构域信息。", self.config.error_log_path)
                    continue
                
                features = data.get('features', [])
                print(f"调试信息: UniProt ID {uniprot_id} 共获取到 {len(features)} 个特征")

                with open(domain_info_file, "a", encoding="utf-8") as domain_f:
                    domain_f.write(f"## {gene} (UniProt ID: {uniprot_id})\n")
                    domain_f.write(f"### 总氨基酸数量: {total_amino_acids}\n\n")
                    
                    found_domains = False
                    for feature in features:
                        feature_type = feature.get('type')
                        if feature_type in ("Domain", "Region"):
                            domain_name = feature.get('description') or feature.get('featureId', '未知结构域')
                            location = feature.get('location', {})
                            begin = location.get('start', {}).get('value')
                            end = location.get('end', {}).get('value')
                            
                            if begin is not None and end is not None:
                                domain_f.write(f"- [{feature_type}] {domain_name}: 序列范围 {begin}-{end}\n")
                                print(f"{gene} 的 {domain_name} 结构域的序列编号范围: {begin}-{end}")
                                found_domains = True
                    if not found_domains:
                        domain_f.write("- 未找到结构域信息。可能原因：\n")
                        domain_f.write("  - 该蛋白质可能没有已知结构域注释\n")
                        domain_f.write("  - UniProt数据库中该条目的注释信息不完整\n")
                        domain_f.write("  - 请检查UniProt ID是否正确或尝试更新UniProt数据\n")
                print(f"已将 {gene} 的结构特征信息写入 {domain_info_file}")

            except requests.exceptions.RequestException as e:
                self.logger.log_error(f"获取 {gene} 的结构域信息失败: {e}", self.config.error_log_path)
            except Exception as e:
                self.logger.log_error(f"处理基因 {gene} 的结构域信息时发生未知错误: {e}", self.config.error_log_path)


    @staticmethod
    def extract_subsequence(fasta_file: Union[str, Path], start: int, end: int) -> Optional[str]:
        """提取序列的指定区域

        Args:
            fasta_file: FASTA文件路径
            start: 起始位置 (1-based)
            end: 结束位置 (1-based)

        Returns:
            提取的子序列或None（如果提取位置无效）
        """
        try:
            _, sequence = CommonUtils.read_fasta(fasta_file, return_header=True)
            if start < 1 or end > len(sequence) or start > end:
                Logger.log_error(f"提取位置超出序列范围或无效 (Start: {start}, End: {end}, Seq Length: {len(sequence)})", ConfigManager().error_log_path)
                return None
            return sequence[start - 1:end]
        except Exception as e:
            Logger.log_error(f"提取子序列失败: {e}", ConfigManager().error_log_path)
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
        sequence_list = list(str(record.seq))
        for pos, aa in zip(mutation_positions, new_amino_acids):
            if 1 <= pos <= len(sequence_list):
                sequence_list[pos - 1] = aa
            else:
                Logger.log_error(f"突变位置 {pos} 超出序列范围。跳过此突变。", ConfigManager().error_log_path)
        mutated_sequence = "".join(sequence_list)
        return SeqRecord(Seq(mutated_sequence), id=record.id, description=record.description)

    @staticmethod
    def read_align_fasta(file_path: Union[str, Path]) -> str:
        """读取比对用的FASTA文件

        Args:
            file_path: 文件路径

        Returns:
            序列字符串
        """
        return str(SeqIO.read(str(file_path), "fasta").seq)

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
        alignments = aligner.align(seq1, seq2)
        if not alignments:
            Logger.log_error("无法找到序列比对。", ConfigManager().error_log_path)
            return 0.0, 0.0

        best_alignment = next(alignments)
        score = best_alignment.score
        
        max_len = max(len(seq1), len(seq2))
        homology = (score / max_len) * 100 if max_len > 0 else 0.0
        
        return score, homology

    def compare_sequences(self, file_paths: List[Union[str, Path]]) -> None:
        """比较多个FASTA文件中的序列

        Args:
            file_paths: FASTA文件路径列表
        """
        if len(file_paths) < 2:
            self.logger.log_error("至少需要两个FASTA文件才能进行比对。", self.config.error_log_path)
            return

        sequences: List[str] = []
        for fp in tqdm(file_paths, desc="Reading alignment files"):
            try:
                sequences.append(self.read_align_fasta(fp))
            except Exception as e:
                self.logger.log_error(f"读取比对文件 {fp} 失败: {e}", self.config.error_log_path)
                return

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
            genes = args.genes
            if not genes:
                self.logger.log_error("未提供有效的基因名称。", self.config.error_log_path)
                return
            
            # 设置基于result/的存储路径
            if args.output_dir is None:
                output_dir, _ = CommonUtils.get_output_dir("result", "protein_sequences", "protein")
            else:
                output_dir = args.output_dir
            
            gene_uniprot_pairs = list(self.fetch_protein_sequences(genes, output_dir))
            if gene_uniprot_pairs:
                self.fetch_domain_information(gene_uniprot_pairs, output_dir)
            else:
                self.logger.log_error("未成功获取任何基因的UniProt ID，无法获取结构域信息。", self.config.error_log_path)

        elif args.command == "extract":
            fasta_file_path = Path(args.fasta_file)
            if not fasta_file_path.is_file():
                self.logger.log_error(f"错误：文件 {fasta_file_path} 不存在。", self.config.error_log_path)
                return
            
            extracted_seq = self.extract_subsequence(fasta_file_path, args.start, args.end)
            if extracted_seq is None:
                return

            header, _ = CommonUtils.read_fasta(fasta_file_path, return_header=True)
            if not isinstance(header, str):
                self.logger.log_error(f"无法从文件 {fasta_file_path} 获取头部信息。", self.config.error_log_path)
                return

            # 设置输出目录
            if args.output_dir is None:
                output_dir = fasta_file_path.parent
            else:
                output_dir = Path(args.output_dir)
                output_dir.mkdir(parents=True, exist_ok=True)

            output_file = CommonUtils.generate_output_filename(fasta_file_path.stem, ".fasta", args.start, args.end)
            output_path = output_dir / output_file
            CommonUtils.save_sequence(header, extracted_seq, output_path)
            print(f"已保存截取序列至 {output_path}")

        elif args.command == "mutate":
            fasta_file_path = Path(args.fasta_file)
            if not fasta_file_path.is_file():
                self.logger.log_error(f"错误：文件 {fasta_file_path} 不存在。", self.config.error_log_path)
                return

            self.validate_input(args.pos, args.aa)
            try:
                record = next(SeqIO.parse(str(fasta_file_path), "fasta"))
            except Exception as e:
                self.logger.log_error(f"解析FASTA文件 {fasta_file_path} 失败: {e}", self.config.error_log_path)
                return

            # 设置输出目录
            if args.output_dir is None:
                output_dir = fasta_file_path.parent
            else:
                output_dir = Path(args.output_dir)
                output_dir.mkdir(parents=True, exist_ok=True)

            mutated_record = self.perform_mutations(record, args.pos, args.aa)
            mutation_info = "".join([f"{p}{a}" for p, a in zip(args.pos, args.aa)])
            output_file = CommonUtils.generate_output_filename(fasta_file_path.stem, ".fasta", mutation_info)
            output_path = output_dir / output_file
            try:
                SeqIO.write(mutated_record, output_path, "fasta")
                print(f"已保存突变序列至 {output_path}")
            except IOError as e:
                self.logger.log_error(f"保存突变序列到文件失败: {output_path} - {e}", self.config.error_log_path)

        elif args.command == "align":
            self.compare_sequences(args.fasta_files)


class PDBProcessor:
    """PDB处理类，负责PDB文件相关操作"""
    def __init__(self, config: ConfigManager, logger: Logger, uniprot_api: 'UniProtAPI'):
        self.config = config
        self.logger = logger
        self.uniprot_api = uniprot_api
        self.exclude_residues = self.parse_exclude_residues()
        self.rdkit_available = 'rdkit' in sys.modules

    def parse_exclude_residues(self) -> List[str]:
        """解析排除残基的配置文件 `exclude_residues.ini`。
        期望配置文件中有一个或多个 [residues] 或类似 section，
        每个 key 对应一个或多个逗号分隔的残基名称。
        """
        config = configparser.ConfigParser()
        config_path = 'exclude_residues.ini'
        if not Path(config_path).exists():
            self.logger.log_error(f"配置文件 {config_path} 未找到。将不排除任何残基。", self.config.error_log_path)
            return []
        
        try:
            config.read(config_path)
        except Exception as e:
            self.logger.log_error(f"读取 {config_path} 配置文件失败: {e}", self.config.error_log_path)
            return []

        exclude_residues = []
        for section in config.sections():
            for key in config[section]:
                exclude_residues.extend([res.strip().upper() for res in config[section][key].split(',') if res.strip()])
        return list(set(exclude_residues))

    @lru_cache(maxsize=None)
    def get_pdb_ids_from_uniprot(self, uniprot_id: str) -> List[str]:
        """通过UniProt ID获取对应的PDB ID列表（带缓存优化）
        使用 UniProt REST API 获取交叉引用信息。

        Args:
            uniprot_id: UniProt ID

        Returns:
            PDB ID列表
        """
        try:
            data = self.uniprot_api.get_uniprot_data(uniprot_id)
            if not data:
                self.logger.log_error(f"无法从UniProt API获取UniProt ID {uniprot_id} 的数据。", self.config.error_log_path)
                return []
            
            pdb_ids = []
            cross_references = data.get('uniProtKBCrossReferences', [])
            for xref in cross_references:
                if xref.get('database') == 'PDB':
                    pdb_ids.append(xref.get('id'))
            return list(set(pdb_ids)) # Return unique PDB IDs
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
        failed_files: List[str] = []

        def download_single_file(pdb_id: str) -> None:
            try:
                local_path = pdbl.retrieve_pdb_file(pdb_id, pdir=str(output_dir_path), file_format='mmCif')
                if local_path and Path(local_path).exists():
                    pass
                else:
                    self.logger.log_error(f"下载CIF文件 {pdb_id} 失败，未返回有效路径或文件不存在。", self.config.error_log_path)
                    failed_files.append(pdb_id)
            except Exception as e:
                error_msg = f"下载CIF文件 {pdb_id} 时出错: {e}"
                self.logger.log_error(error_msg, self.config.error_log_path)
                failed_files.append(pdb_id)

        CommonUtils.parallel_executor(download_single_file, pdb_ids, description="Downloading PDB files")

        if failed_files:
            failed_log_path = output_dir_path / 'download_failed.txt'
            with open(failed_log_path, 'a', encoding='utf-8') as f:
                for failed_file_id in failed_files:
                    f.write(f"{failed_file_id}.cif\n")
            print(f"部分CIF文件下载失败，详情请查看 {failed_log_path}")

    def extract_ligands_from_pdb(self, folder_path: Union[str, Path]) -> Dict[str, List[str]]:
        """从.cif文件中提取配体信息

        Args:
            folder_path: 文件夹路径

        Returns:
            配体到PDB ID的映射字典
        """
        parser = MMCIFParser()
        ligand_pdb_dict: Dict[str, List[str]] = {}

        folder_path = Path(folder_path)
        if not folder_path.is_dir():
            self.logger.log_error(f"指定的PDB文件夹路径不存在或不是目录: {folder_path}", self.config.error_log_path)
            return {}

        cif_files = list(folder_path.glob('*.cif'))
        if not cif_files:
            self.logger.log_error(f"在 {folder_path} 中未找到任何.cif文件进行配体提取。", self.config.error_log_path)
            return {}

        for file_path in tqdm(cif_files, desc="Extracting ligands"):
            pdb_id = file_path.stem
            try:
                structure = parser.get_structure(pdb_id, str(file_path))
                ligands_in_pdb: List[str] = []
                for model in structure:
                    for chain in model:
                        for residue in chain:
                            if residue.id[0].startswith('H_') and residue.resname.strip().upper() not in self.exclude_residues:
                                ligands_in_pdb.append(residue.resname.strip().upper())
                
                for ligand in set(ligands_in_pdb):
                    ligand_pdb_dict.setdefault(ligand, []).append(pdb_id)

            except Exception as e:
                error_msg = f"解析 {file_path.name} 时出错: {e}"
                self.logger.log_error(error_msg, self.config.error_log_path)
                continue
        return ligand_pdb_dict

    @staticmethod
    def write_ligand_info_to_md(ligand_pdb_dict: Dict[str, List[str]], output_dir: Union[str, Path]) -> None:
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
                for ligand, pdb_ids in sorted(ligand_pdb_dict.items()):
                    pdb_id_str = ', '.join(sorted(list(set(pdb_ids))))
                    outfile.write(f'| {ligand} | {pdb_id_str} |\n')
            print(f'结果已写入 {output_file_path}')
        except IOError as e:
            Logger.log_error(f"写入配体信息到MD文件失败: {output_file_path} - {e}", ConfigManager().error_log_path)
            raise

    def move_files_based_on_ligands(self, folder_path: Union[str, Path], ligand_pdb_dict: Dict[str, List[str]]) -> None:
        """根据配体信息移动CIF文件到'no_ligand'或'with_ligands'子文件夹

        Args:
            folder_path: 包含CIF文件的文件夹路径
            ligand_pdb_dict: 配体到PDB ID的映射字典
        """
        folder_path = Path(folder_path)
        no_ligand_dir = folder_path / 'no_ligand'
        with_ligands_dir = folder_path / 'with_ligands'

        no_ligand_dir.mkdir(exist_ok=True)
        with_ligands_dir.mkdir(exist_ok=True)

        pdb_ids_with_ligands = set()
        for pdb_ids in ligand_pdb_dict.values():
            pdb_ids_with_ligands.update(pdb_ids)

        cif_files = list(folder_path.glob('*.cif'))
        for src_path in tqdm(cif_files, desc="Moving PDB files"):
            pdb_id = src_path.stem
            
            if pdb_id in pdb_ids_with_ligands:
                dst_path = with_ligands_dir / src_path.name
                target_dir_name = 'with_ligands'
            else:
                dst_path = no_ligand_dir / src_path.name
                target_dir_name = 'no_ligand'
            
            try:
                if src_path != dst_path:
                    shutil.move(src_path, dst_path)
            except Exception as e:
                self.logger.log_error(f"移动文件 {src_path.name} 到 {target_dir_name} 时出错: {e}", self.config.error_log_path)

    def download_ligand_json(self, unique_ligands: Set[str], output_dir: Union[str, Path]) -> None:
        """通过API并行查询配体的JSON文件并下载保存到json子文件夹

        Args:
            unique_ligands: 唯一配体集合
            output_dir: 输出目录
        """
        json_dir = Path(output_dir) / 'json'
        json_dir.mkdir(exist_ok=True)

        def download_single_ligand_json(ligand: str) -> None:
            session = self.uniprot_api.session
            try:
                url = f"https://data.rcsb.org/rest/v1/core/chemcomp/{ligand}"
                response = session.get(url, timeout=10)
                response.raise_for_status()
                
                data = response.json()
                json_path = json_dir / f'{ligand}.json'
                with open(json_path, 'w', encoding='utf-8') as json_file:
                    json.dump(data, json_file, ensure_ascii=False, indent=4)
            except requests.exceptions.HTTPError as e:
                if e.response.status_code == 404:
                    self.logger.log_error(f"配体 {ligand} 未在RCSB PDB找到 (404 Not Found)。", self.config.error_log_path)
                else:
                    self.logger.log_error(f"请求配体 {ligand} 失败，状态码: {e.response.status_code} - {e.response.text}", self.config.error_log_path)
            except requests.exceptions.RequestException as e:
                self.logger.log_error(f"下载配体 {ligand} 时发生网络错误: {e}", self.config.error_log_path)
            except json.JSONDecodeError:
                self.logger.log_error(f"下载配体 {ligand} 时API返回无效JSON。", self.config.error_log_path)
            except Exception as e:
                self.logger.log_error(f"下载配体 {ligand} 时出错: {e}", self.config.error_log_path)

        CommonUtils.parallel_executor(download_single_ligand_json, list(unique_ligands), description="Downloading ligand JSONs")

    def write_chemical_info_to_md(self, unique_ligands: Set[str], output_dir: Union[str, Path]) -> None:
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
                for ligand in tqdm(sorted(list(unique_ligands)), desc="Writing chemical info"):
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
                    except FileNotFoundError:
                        self.logger.log_error(f"JSON file for {ligand} not found in {json_dir}", self.config.error_log_path)
                    except Exception as e:
                        self.logger.log_error(f"读取或解析配体 {ligand} 的JSON文件失败: {e}", self.config.error_log_path)
            print(f"MD file {output_md_file} has been created successfully.")
        except IOError as e:
            self.logger.log_error(f"创建或写入 chemical_components_info.md 文件失败: {e}", self.config.error_log_path)
            raise

    def calculate_similarity_and_rank_pdbs(self, user_smiles: str, output_dir: str, top_n: int = 5) -> List[Tuple[str, float]]:
        """计算用户输入的小分子SMILES与PDB结构中配体的结构相似性，并排序PDB结构

        Args:
            user_smiles: 用户输入的小分子SMILES字符串
            output_dir: 输出目录路径
            top_n: 返回前N个最高相似性的PDB结构

        Returns:
            排序后的PDB ID和相似性分数的元组列表
        """
        if not self.rdkit_available:
            self.logger.log_error("RDKit库未正确安装，无法计算结构相似性。", self.config.error_log_path)
            return []

        # 读取chemical_components_info.md文件获取配体SMILES信息
        md_file_path = Path(output_dir) / 'chemical_components_info.md'
        if not md_file_path.exists():
            self.logger.log_error(f"配体信息文件 {md_file_path} 不存在，无法获取配体SMILES信息。", self.config.error_log_path)
            return []

        try:
            # 读取用户输入的分子
            user_mol = Chem.MolFromSmiles(user_smiles)
            if user_mol is None:
                self.logger.log_error(f"无效的SMILES字符串: {user_smiles}", self.config.error_log_path)
                return []
            user_fp = AllChem.GetMorganFingerprintAsBitVect(user_mol, 2, nBits=2048)

            # 读取配体信息
            ligand_smiles_dict = {}
            with open(md_file_path, 'r', encoding='utf-8') as f:
                lines = f.readlines()
            for line in lines[2:]:  # 跳过表头和分隔线
                parts = [p.strip() for p in line.strip().split('|') if p.strip()]
                if len(parts) >= 5 and parts[4]:  # 检查是否有SMILES信息
                    ligand_id = parts[0]
                    smiles = parts[4]
                    ligand_smiles_dict[ligand_id] = smiles

            # 读取pdb_ligand.md获取配体-PDB映射
            pdb_ligand_file = Path(output_dir) / 'pdb_ligand.md'
            ligand_pdb_mapping = {}
            if pdb_ligand_file.exists():
                with open(pdb_ligand_file, 'r', encoding='utf-8') as f:
                    lines = f.readlines()
                for line in lines[2:]:  # 跳过表头和分隔线
                    parts = [p.strip() for p in line.strip().split('|') if p.strip()]
                    if len(parts) >= 2:
                        ligand_id = parts[0]
                        pdb_ids = [pdb_id.strip() for pdb_id in parts[1].split(',')]
                        ligand_pdb_mapping[ligand_id] = pdb_ids

            # 计算相似性并排序
            similarity_results = []
            processed_pdbs = set()  # 避免重复的PDB ID

            for ligand_id, pdb_ids in ligand_pdb_mapping.items():
                if ligand_id in ligand_smiles_dict:
                    ligand_smiles = ligand_smiles_dict[ligand_id]
                    try:
                        ligand_mol = Chem.MolFromSmiles(ligand_smiles)
                        if ligand_mol is not None:
                            ligand_fp = AllChem.GetMorganFingerprintAsBitVect(ligand_mol, 2, nBits=2048)
                            similarity = DataStructs.TanimotoSimilarity(user_fp, ligand_fp)
                            
                            for pdb_id in pdb_ids:
                                if pdb_id not in processed_pdbs:
                                    similarity_results.append((pdb_id, similarity))
                                    processed_pdbs.add(pdb_id)
                    except Exception as e:
                        self.logger.log_error(f"计算配体 {ligand_id} 相似性时出错: {e}", self.config.error_log_path)
                        continue

            # 根据相似性分数排序，取前N个
            similarity_results.sort(key=lambda x: x[1], reverse=True)
            top_results = similarity_results[:top_n]
            
            # 将结果写入文件
            result_file = Path(output_dir) / 'structure_similarity_ranking.md'
            with open(result_file, 'w', encoding='utf-8') as f:
                f.write("# 分子结构相似性排序结果\n\n")
                f.write(f"用户输入的SMILES: {user_smiles}\n\n")
                f.write("| 排名 | PDB ID | 结构相似性分数 (Tanimoto系数) |\n")
                f.write("| --- | --- | --- |\n")
                for i, (pdb_id, score) in enumerate(top_results, 1):
                    f.write(f"| {i} | {pdb_id} | {score:.4f} |\n")
            
            print(f"已将结构相似性排序结果写入 {result_file}")
            print(f"前{min(len(top_results), top_n)}个最相似的PDB结构：")
            for i, (pdb_id, score) in enumerate(top_results, 1):
                print(f"{i}. {pdb_id}: {score:.4f}")
            
            return top_results
        except Exception as e:
            self.logger.log_error(f"计算结构相似性时出错: {e}", self.config.error_log_path)
            return []

    def process(self, uniprot_id: str, output_dir: str, user_smiles: str = None) -> None:
        """处理PDB文件下载和分析

        Args:
            uniprot_id: UniProt蛋白质编号
            output_dir: 输出目录路径
            user_smiles: 可选，用户输入的小分子SMILES字符串，用于计算结构相似性
        """
        output_dir_path = Path(output_dir)
        output_dir_path.mkdir(parents=True, exist_ok=True)

        if not uniprot_id:
            self.logger.log_error("缺少 'uniprot' ID。", self.config.error_log_path)
            return
        
        print(f"\n--- 开始处理 PDB 任务 (UniProt ID: {uniprot_id}) ---")
        
        pdb_ids = self.get_pdb_ids_from_uniprot(uniprot_id)
        if pdb_ids:
            print(f"为 UniProt ID {uniprot_id} 找到 PDB IDs: {', '.join(pdb_ids)}")
            self.download_pdb_files(pdb_ids, output_dir)
        else:
            self.logger.log_error(f"未找到 UniProt ID {uniprot_id} 对应的PDB ID。", self.config.error_log_path)
            return

        downloaded_cif_files = list(Path(output_dir).glob('*.cif'))
        if not downloaded_cif_files:
            self.logger.log_error(f"在 {output_dir} 中未找到下载的CIF文件，跳过后续配体处理。", self.config.error_log_path)
            return

        ligand_pdb_dict = self.extract_ligands_from_pdb(output_dir)
        
        if ligand_pdb_dict:
            self.write_ligand_info_to_md(ligand_pdb_dict, output_dir)
        else:
            print(f"在 {output_dir} 中未找到任何配体信息。")
        
        self.move_files_based_on_ligands(output_dir, ligand_pdb_dict)
        md_file_path = Path(output_dir) / 'pdb_ligand.md'
        unique_ligands: Set[str] = set()
        if md_file_path.exists():
            try:
                with open(md_file_path, 'r', encoding='utf-8') as f:
                    lines = f.readlines()
                for line in lines[2:]:
                    parts = [p.strip() for p in line.strip().split('|') if p.strip()]
                    if len(parts) >= 1:
                        unique_ligands.add(parts[0])
            except IOError as e:
                self.logger.log_error(f"读取 {md_file_path} 失败: {e}", self.config.error_log_path)
        else:
            self.logger.log_error(f"配体信息文件 {md_file_path} 不存在，无法获取唯一配体列表。", self.config.error_log_path)

        if unique_ligands:
            print(f"找到唯一配体: {', '.join(sorted(list(unique_ligands)))}")
            self.download_ligand_json(unique_ligands, output_dir)
            self.write_chemical_info_to_md(unique_ligands, output_dir)
        else:
            print("未找到任何配体，跳过配体JSON下载和坐标提取。")

        # 如果提供了用户SMILES，计算结构相似性并排序PDB结构
        if user_smiles:
            self.logger.log_info(f"开始计算分子结构相似性，用户输入的SMILES: {user_smiles}", self.config.info_log_path)
            self.calculate_similarity_and_rank_pdbs(user_smiles, str(output_dir))

        print('\nPDB任务处理完成。')


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
            'User-Agent': 'PyPDA/0.1.0 (https://gitee.com/coding_playground/py-pda; dengzho5068@foxmail.com)'
        })
        return session

    def search_uniprot_by_name(self, name: str, organism_id: int = 9606) -> Optional[str]:
        """
        根据蛋白质名称或基因名称搜索UniProt，并限定种属。
        Args:
            name: 蛋白质或基因名称。
            organism_id: 限制的物种ID (默认9606为人类)。
        Returns:
            找到的第一个UniProt ID，或None。
        """
        params = {
            "query": f"gene_exact:{name} AND organism_id:{organism_id}",
            "format": "json",
            "fields": "accession"
        }
        search_url = f"{self.api_base_url}search"
        try:
            print(f"正在搜索 UniProt ID for '{name}' (Taxon ID: {organism_id})...")
            response = self.session.get(search_url, params=params, timeout=30)
            response.raise_for_status()
            data = response.json()
            
            if data and 'results' in data and data['results']:
                accession = data['results'][0].get('primaryAccession')
                if accession:
                    print(f"为 '{name}' 找到 UniProt ID: {accession}")
                    return accession
            return None
        except requests.exceptions.RequestException as e:
            self.logger.log_error(f"搜索 UniProt ID for '{name}' 失败: {e}", self.config.error_log_path)
            return None
        except json.JSONDecodeError:
            self.logger.log_error(f"UniProt API搜索返回无效JSON for '{name}': {response.text[:200]}...", self.config.error_log_path)
            return None
        except Exception as e:
            self.logger.log_error(f"处理 UniProt 搜索结果失败 for '{name}': {e}", self.config.error_log_path)
            return None

    def get_uniprot_data(self, accession: str) -> Optional[Dict[str, Any]]:
        """根据UniProt编号从API获取数据

        Args:
            accession: UniProt蛋白质编号

        Returns:
            如果成功获取，返回蛋白质信息的字典；否则返回None。
        """
        url = f'{self.api_base_url}{accession}'
        try:
            print(f"正在从 {url} 获取数据...")
            response = self.session.get(url, timeout=30)
            response.raise_for_status()
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
            self.logger.log_error(f"获取UniProt数据失败: {e} (UniProt API - {url})", self.config.error_log_path)
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

        # 4. 功能和活性注释
        comments = data.get('comments', [])
        comment_types = ['FUNCTION', 'CATALYTIC ACTIVITY', 'COFACTOR', 'ACTIVITY REGULATION',
                         'TISSUE SPECIFICITY', 'SUBCELLULAR LOCATION', 'PTM', 'SIMILARITY', 'DISEASE', 'INTERACTION']
        comment_type_map = {ct.lower(): ct for ct in comment_types}
        info['comments'] = {ct: [] for ct in comment_types}
        
        for comment in comments:
            ct_lower = comment.get('commentType', '').lower()
            if ct_lower in comment_type_map:
                ct = comment_type_map[ct_lower]
                info['comments'][ct].append(comment)
                if ct == 'INTERACTION':
                    info['interactions'].extend(comment.get('interactions', []))

        # 5. 蛋白质特征
        features = data.get('features', [])
        info['features'] = {
            'counts': {},
            'detailed': {}
        }
        detailed_feature_types = ["Chain", "Region", "Active site", "Binding site", "Modified residue", "Mutagenesis", "Domain"]
        for dt in detailed_feature_types:
            info['features']['detailed'][dt] = []

        for feature in features:
            ft = feature.get('type')
            info['features']['counts'][ft] = info['features']['counts'].get(ft, 0) + 1
            if ft == "Region" and feature.get('description') == "Domain":
                info['features']['detailed']['Domain'].append(feature)
            elif ft in info['features']['detailed']:
                info['features']['detailed'][ft].append(feature)

        # 6. 蛋白质相互作用(已在4. 功能和活性注释中提取)

        # 7. 关键词
        info['keywords'] = [kw.get('name') for kw in data.get('keywords', []) if kw.get('name')]

        # 8. 参考文献
        info['references'] = data.get('references', [])

        # 9. 交叉引用
        cross_references = data.get('uniProtKBCrossReferences', [])
        important_databases = ["PDB", "DrugBank", "GO", "Reactome", "HGNC", "GeneID", "KEGG", "AlphaFoldDB", "IntAct"]
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
    def generate_md_report(info: ProteinInfo, md_filename: Union[str, Path]) -> None:
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
                                 'DISEASE',]
                
                found_any_comment = False
                for ct in comment_order:
                    if comments_data.get(ct):
                        found_any_comment = True
                        f.write(f"### {ct.replace('_', ' ').title()}\n")
                        for comment in comments_data[ct]:
                            texts = []
                            if 'texts' in comment and isinstance(comment['texts'], list):
                                for text_item in comment['texts']:
                                    if isinstance(text_item, dict) and 'value' in text_item:
                                        texts.append(text_item['value'])
                                    elif isinstance(text_item, str):
                                        texts.append(text_item)
                            elif 'text' in comment:
                                if isinstance(comment['text'], dict) and 'value' in comment['text']:
                                    texts.append(comment['text']['value'])
                                elif isinstance(comment['text'], str):
                                    texts.append(comment['text'])

                            if texts:
                                f.write(f"- {'; '.join(texts)}\n")
                            
                            if ct == 'CATALYTIC ACTIVITY':
                                reaction = comment.get('reaction')
                                if reaction:
                                    f.write(f"  - **Reaction**: {reaction.get('name', 'N/A')}\n")
                            elif ct == 'DISEASE':
                                disease_name = comment.get('disease', {}).get('diseaseName', 'N/A')
                                acronym = comment.get('disease', {}).get('acronym', 'N/A')
                                f.write(f"  - **Disease**: {disease_name} ({acronym})\n")
                            elif ct == 'SUBCELLULAR LOCATION':
                                for location in comment.get('subcellularLocations', []):
                                    f.write(f"  - **Location**: {location.get('location', {}).get('value', 'N/A')}\n")
                                    for topology in location.get('topologies', []):
                                        f.write(f"    - **Topology**: {topology.get('value', 'N/A')}\n")
                                    for orientation in location.get('orientations', []):
                                        f.write(f"    - **Orientation**: {orientation.get('value', 'N/A')}\n")
                            elif ct == 'COFACTOR':
                                cofactors = comment.get('cofactors', [])
                                if cofactors:
                                    f.write("  - **Cofactors**:\n")
                                    for cofactor in cofactors:
                                        cofactor_name = cofactor.get('name', 'N/A')
                                        chebi_id = cofactor.get('cofactorCrossReference', {}).get('id', 'N/A')
                                        if chebi_id != 'N/A':
                                            f.write(f"    - {cofactor_name} (ChEBI ID: {chebi_id})\n")
                                        else:
                                            f.write(f"    - {cofactor_name}\n")
                                else:
                                    f.write("  - 未找到详细辅因子信息。\n")
          
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
                detailed_feature_order = ["Chain", "Region", "Active site", "Binding site", "Modified residue", "Mutagenesis", "Domain"]
                for dt in detailed_feature_order:
                    if features['detailed'].get(dt):
                        found_any_detailed_feature = True
                        f.write(f"#### {dt}\n")
                        for feature in features['detailed'][dt]:
                            desc = feature.get('description', 'N/A')
                            if isinstance(desc, dict):
                                desc = desc.get('value', 'N/A')

                            loc = feature.get('location', {})
                            begin = loc.get('start', {}).get('value', 'N/A')
                            end = loc.get('end', {}).get('value', 'N/A')
                            
                            extra_info = []
                            if dt == 'Mutagenesis':
                                alt_seq_data = feature.get('alternativeSequence', {})
                                original_seq = alt_seq_data.get('originalSequence')
                                alternative_seqs = alt_seq_data.get('alternativeSequences', [])
                                
                                if original_seq:
                                    mutation_str = f"原始: {original_seq}"
                                    if alternative_seqs:
                                        mutation_str += f" -> 突变: {', '.join(alternative_seqs)}"
                                    extra_info.append(mutation_str)
                            
                            feature_line = f"- {desc} (位置: {begin}-{end})"
                            if extra_info:
                                feature_line += f", {', '.join(extra_info)}" 

                            f.write(feature_line + "\n")
                if not found_any_detailed_feature:
                    f.write("- 未找到详细特征信息。\n")
                f.write("\n")

                # 6. 蛋白质相互作用
                f.write("## 6. 蛋白质相互作用\n")
                interactions = info['interactions']
                if interactions:
                    for i, interaction in enumerate(interactions, 1):
                        interactant_one = interaction.get('interactantOne', {})
                        interactant_two = interaction.get('interactantTwo', {})
                        main_protein_id = interactant_one.get('uniProtKBAccession', 'N/A')
                        partner_id = interactant_two.get('uniProtKBAccession', 'N/A')
                        partner_name = interactant_two.get('geneName', interactant_two.get('name', 'N/A'))
                        organism_differ = "是" if interaction.get('organismDiffer', False) else "否"
                        f.write(f"- 相互作用 {i}: {main_protein_id} <-> **{partner_name}** (UniProt ID: {partner_id})\n")
                else:
                    f.write("- 未找到相互作用信息\n")
                f.write("\n")

                # 7. 关键词
                f.write("## 7. 关键词\n")
                keywords = info['keywords']
                if keywords:
                    f.write(f"- {', '.join(sorted(keywords))}\n")
                else:
                    f.write("- 未找到关键词信息\n")
                f.write("\n")

                # 8. 参考文献
                f.write("## 8. 参考文献\n")
                references = info['references']
                if references:
                    for i, ref in enumerate(references[:10], 1):
                        citation = ref.get('citation', {})
                        authors_list = citation.get('authors', [])
                        authors = ", ".join(authors_list) if authors_list else 'N/A'
                        title = citation.get('title', 'N/A')
                        journal = citation.get('journal', 'N/A')
                        year = citation.get('publicationDate', 'N/A')[:4]
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
                for db in sorted(cross_refs.keys()):
                    if cross_refs[db]:
                        found_any_xref = True
                        f.write(f"### {db}\n")
                        ids = [xref.get('id') for xref in cross_refs[db] if xref.get('id')]
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
                    formatted_seq = '\n'.join([seq[i:i+80] for i in range(0, len(seq), 80)])
                    f.write("- **氨基酸序列**:\n\n```\n{}\n```\n".format(formatted_seq))
                else:
                    f.write("- 未找到氨基酸序列信息。\n")

            print(f"Markdown报告已生成至 {md_filename}")
        except IOError as e:
            Logger.log_error(f"生成Markdown报告失败: {md_filename} - {e}", ConfigManager().error_log_path)
            raise


class PypdaApp:
    """Pypda应用主类，负责命令行参数解析和工具调度"""
    def __init__(self):
        self.config = ConfigManager()
        self.logger = Logger()
        self.uniprot_api = UniProtAPI(self.config, self.logger)
        self.seq_processor = SequenceProcessor(self.config, self.logger, self.uniprot_api)
        self.pdb_processor = PDBProcessor(self.config, self.logger, self.uniprot_api)

    def _get_human_uniprot_id_and_handle_error(self, protein_name: str, task_context: str) -> Optional[str]:
        """
        尝试从蛋白质名称获取人类UniProt ID。
        如果失败，则记录特定于任务的错误并返回None，表示调用者应停止当前任务。
        Args:
            protein_name: 蛋白质或基因名称。
            task_context: 当前任务的描述（例如“PDB处理”或“UniProt数据获取”）。
        Returns:
            成功解析的UniProt ID，否则为None。
        """
        uniprot_id = self.uniprot_api.search_uniprot_by_name(protein_name, organism_id=9606)
        if not uniprot_id:
            self.logger.log_error(f"无法为 '{protein_name}' (人类) 获取UniProt ID，跳过{task_context}。", self.config.error_log_path)
        return uniprot_id

    def _setup_seq_parser(self, subparsers: argparse._SubParsersAction) -> None:
        """设置序列分析工具的子命令解析器"""
        seq_parser = subparsers.add_parser(
            "seq", 
            help="蛋白质序列处理工具，支持序列获取、提取、突变和比对。",
            formatter_class=argparse.RawTextHelpFormatter
        )
        seq_subparsers = seq_parser.add_subparsers(dest="command", required=True, help="序列处理命令")

        # seq fetch命令
        fetch_parser = seq_subparsers.add_parser(
            "fetch", 
            help="批量获取蛋白质序列和结构域信息。",
            description="""
            从UniProt获取指定基因的人类蛋白质全长序列，并提取结构域信息。
            基因名称用空格分隔输入（例如：BRCA1 TP53 EGFR）。
            """
        )
        fetch_parser.add_argument("genes", nargs='+', type=str, help="要下载的基因名称，用空格分隔（例如：BRCA1 TP53 EGFR）。")
        fetch_parser.add_argument(
            "output_dir", 
            type=str,
            default=None, 
            nargs='?',
            help="输出目录，用于保存FASTA序列文件和结构域信息报告 (默认: result/protein_sequences/基因名_YYYYMMDD_HHMMSS)。"
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
        extract_parser.add_argument(
            "--output_dir", 
            type=str, 
            default=None, 
            help="输出目录，用于保存提取的子序列文件 (默认: 输入文件所在目录)。"
        )

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
        mut_parser.add_argument(
            "--output_dir", 
            type=str, 
            default=None, 
            help="输出目录，用于保存突变后的序列文件 (默认: 输入文件所在目录)。"
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
            根据提供的蛋白质名称（或基因名称），搜索对应的人类蛋白质UniProt ID，
            然后自动下载相关PDB文件 (mmCIF格式)，提取蛋白质中的配体信息，
            生成报告，并根据是否含有配体将PDB文件分类。
            同时，下载配体化学信息并提取配体坐标。
            使用--smile参数可根据输入的小分子SMILES计算结构相似性并排序PDB结构。
            """)
        pdb_parser.add_argument("protein_name", type=str, help="蛋白质名称或基因名称，例如: BRCA1。")
        pdb_parser.add_argument(
            "output_dir", 
            type=str, 
            default=None, 
            help="输出目录，用于保存PDB文件和分析结果 (默认: result/pdb_output/蛋白质名_YYYYMMDD_HHMMSS)。")
        pdb_parser.add_argument(
            "--smile",
            type=str,
            default=None,
            help="输入小分子的SMILES字符串，用于计算与PDB结构中配体的结构相似性并排序。")

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
            根据蛋白质名称（或基因名称），搜索对应的人类蛋白质UniProt ID，
            从UniProt REST API获取完整的蛋白质信息，
            并将其保存为JSON文件和Markdown格式的分析报告。
            """
        )
        uniprot_fetch_parser.add_argument("protein_name", type=str, help="蛋白质名称或基因名称，例如: TP53。")
        uniprot_fetch_parser.add_argument(
            "-o", "--output_dir", 
            type=str, 
            default=None, 
            help="保存JSON和Markdown报告的输出目录 (默认: result/uniprot_reports/蛋白质名_YYYYMMDD_HHMMSS)。"
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
            formatter_class=argparse.RawTextHelpFormatter
        )
        
        parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1.0')

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
                uniprot_id = self._get_human_uniprot_id_and_handle_error(args.protein_name, "PDB处理")
                if not uniprot_id:
                    return

                # 设置基于result/的存储路径
                if args.output_dir is None:
                    output_dir, _ = CommonUtils.get_output_dir("result", "pdb_output", args.protein_name)
                else:
                    output_dir = args.output_dir

                self.pdb_processor.process(uniprot_id=uniprot_id, output_dir=output_dir, user_smiles=args.smile)
            elif args.tool == "uniprot":
                if args.command == "fetch":
                    # 设置基于result/的存储路径
                    if args.output_dir is None:
                        output_dir_path_str, _ = CommonUtils.get_output_dir("result", "uniprot_reports", args.protein_name)
                        output_dir_path = Path(output_dir_path_str)
                    else:
                        output_dir_path = Path(args.output_dir)
                    output_dir_path.mkdir(parents=True, exist_ok=True)

                    uniprot_id = self._get_human_uniprot_id_and_handle_error(args.protein_name, "UniProt数据获取")
                    if not uniprot_id:
                        return

                    data = self.uniprot_api.get_uniprot_data(uniprot_id)
                    if not data:
                        self.logger.log_error("无法获取UniProt数据。", self.config.error_log_path)
                        return
                    json_filename = CommonUtils.save_to_json(data, uniprot_id, str(output_dir_path))
                    analyzer = ProteinAnalyzer()
                    protein_info = analyzer.extract_protein_info(data)
                    md_filename = str(output_dir_path / (Path(json_filename).stem + '.md'))
                    ReportGenerator.generate_md_report(protein_info, md_filename)
                elif args.command == "analyze":
                    data = CommonUtils.load_from_json(args.file)
                    analyzer = ProteinAnalyzer()
                    protein_info = analyzer.extract_protein_info(data)
                    md_filename = str(Path(args.file).with_suffix('.md'))
                    ReportGenerator.generate_md_report(protein_info, md_filename)
        except Exception as e:
            self.logger.log_error(f"应用程序运行过程中发生未捕获的错误: {e}", self.config.error_log_path)

if __name__ == "__main__":
    app = PypdaApp()
    app.run()