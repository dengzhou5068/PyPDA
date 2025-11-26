#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
序列处理模块
"""
import argparse
import io
import os
from pathlib import Path
from typing import List, Tuple, Optional, Any, Union

from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm

from config.config_manager import ConfigManager
from logger.logger import Logger
from utils.common_utils import CommonUtils


class SequenceProcessor:
    """序列处理类，负责蛋白质序列相关操作"""
    def __init__(self, config: ConfigManager, logger: Logger, uniprot_api: Any):
        self.config = config
        self.logger = logger
        self.uniprot_api = uniprot_api
    def fetch_protein_sequences(self, genes: List[str], output_dir: str) -> List[Tuple[str, str]]:
        """从UniProt获取蛋白质序列

        Args:
            genes: 基因名称列表
            output_dir: 输出目录

        Yields:
            基因名称和对应的UniProt ID元组
        """
        output_dir_path = Path(output_dir)
        output_dir_path.mkdir(parents=True, exist_ok=True)

        gene_uniprot_pairs = []
        for gene in tqdm(genes, desc="Fetching sequences"):
            # 定义查询策略列表，从最精确到最宽松
            query_strategies = [
                f"gene_exact:{gene} AND organism_id:9606",  # 精确匹配基因名
                f"gene:{gene} AND organism_id:9606",         # 宽松匹配基因名
                f"{gene} AND organism_id:9606"               # 一般搜索，匹配任何字段
            ]
            
            found = False
            for i, query in enumerate(query_strategies):
                params = {
                    "query": query,
                    "format": "fasta",
                    "fields": "accession,sequence"
                }
                try:
                    print(f"尝试查询策略 {i+1} 用于基因 {gene}: {query}")
                    response = self.uniprot_api.session.get(
                        f"{self.config.uniprot_api_base_url}search",
                        params=params,
                        timeout=30
                    )
                    response.raise_for_status()

                    if response.text.strip():
                        output_file = output_dir_path / f"{gene}.fasta"
                        try:
                            first_seq_record = next(SeqIO.parse(io.StringIO(response.text.strip()), "fasta"))
                            SeqIO.write(first_seq_record, output_file, "fasta")
                            
                            print(f"已将 {gene} 的第一个全长蛋白质序列保存到 {output_file}")
                            uniprot_id = first_seq_record.id.split('|')[1] if '|' in first_seq_record.id else first_seq_record.id
                            
                            # 验证获取到的序列是否真正匹配
                            seq_desc = first_seq_record.description
                            if gene.lower() in seq_desc.lower() or uniprot_id:
                                found = True
                                gene_uniprot_pairs.append((gene, uniprot_id))
                                break
                            else:
                                print(f"警告：获取到的序列描述似乎不匹配 {gene}，尝试下一个查询策略")
                                continue
                        except StopIteration:
                            print(f"警告：UniProt API返回无效FASTA格式，尝试下一个查询策略")
                            continue
                except Exception as e:
                    print(f"警告：查询策略 {i+1} 失败: {e}，尝试下一个查询策略")
                    continue
            
            if not found:
                self.logger.log_error(f"未找到基因 {gene} 对应的蛋白质序列", self.config.error_log_path)
                print(f"未找到基因 {gene} 的蛋白质序列，已尝试所有查询策略")
                
                # 尝试通过search_uniprot_by_name方法获取UniProt ID
                try:
                    print(f"尝试使用search_uniprot_by_name方法查询 {gene}")
                    uniprot_id = self.uniprot_api.search_uniprot_by_name(gene, 9606)
                    if uniprot_id:
                        # 如果找到UniProt ID，则直接获取该ID的序列
                        params = {
                            "query": f"accession:{uniprot_id}",
                            "format": "fasta",
                            "fields": "accession,sequence"
                        }
                        response = self.uniprot_api.session.get(
                            f"{self.config.uniprot_api_base_url}search",
                            params=params,
                            timeout=30
                        )
                        response.raise_for_status()
                        
                        if response.text.strip():
                            output_file = output_dir_path / f"{gene}.fasta"
                            first_seq_record = next(SeqIO.parse(io.StringIO(response.text.strip()), "fasta"))
                            SeqIO.write(first_seq_record, output_file, "fasta")
                            print(f"通过UniProt ID {uniprot_id} 成功获取 {gene} 的序列")
                            found = True
                            gene_uniprot_pairs.append((gene, uniprot_id))
                except Exception as e:
                    print(f"通过search_uniprot_by_name方法查询失败: {e}")
                    pass
        
        return gene_uniprot_pairs

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

            except Exception as e:
                self.logger.log_error(f"获取 {gene} 的结构域信息失败: {e}", self.config.error_log_path)

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
                Logger.log_error(f"提取位置超出序列范围或无效 (Start: {start}, End: {end}, Seq Length: {len(sequence)})")
                return None
            return sequence[start - 1:end]
        except Exception as e:
            Logger.log_error(f"提取子序列失败: {e}")
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
            Logger.log_error("突变位置和新氨基酸的数量必须相同。")
            return False
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
                Logger.log_error(f"突变位置 {pos} 超出序列范围。跳过此突变。")
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
        try:
            aligner = PairwiseAligner()
            # 设置为全局比对模式并限制只返回最佳对齐
            aligner.mode = 'global'
            
            # 直接使用score方法获取最佳得分，避免生成所有可能的对齐
            score = aligner.score(seq1, seq2)
            
            # 计算同源性百分比
            max_len = max(len(seq1), len(seq2))
            homology = (score / max_len) * 100 if max_len > 0 else 0.0
            
            return score, homology
        except Exception as e:
            Logger.log_error(f"序列比对时发生错误: {str(e)}")
            return 0.0, 0.0

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
            
            gene_uniprot_pairs = self.fetch_protein_sequences(genes, output_dir)
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

            if not self.validate_input(args.pos, args.aa):
                return
            
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
