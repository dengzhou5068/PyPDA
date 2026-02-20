#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PyPDA 命令行入口脚本
"""

import sys
import os
import argparse
from pathlib import Path

# 将当前目录添加到Python路径
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from config.config_manager import ConfigManager
from logger.logger import Logger
from sequence.sequence_processor import SequenceProcessor
from pdb.pdb_processor import PDBProcessor
from uniprot.uniprot_api import UniProtAPI
from uniprot.protein_analyzer import ProteinAnalyzer
from report.report_generator import ReportGenerator
from utils.common_utils import CommonUtils


class PypdaApp:
    """Pypda应用主类，负责命令行参数解析和工具调度"""
    def __init__(self):
        self.config = ConfigManager()
        self.logger = Logger()
        self.uniprot_api = UniProtAPI(self.config, self.logger)
        self.seq_processor = SequenceProcessor(self.config, self.logger, self.uniprot_api)
        self.pdb_processor = PDBProcessor(self.config, self.logger, self.uniprot_api)

    def _get_human_uniprot_id_and_handle_error(self, protein_name: str, task_context: str) -> str:
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
            "--output_dir", 
            type=str,
            default=None, 
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
            formatter_class=argparse.RawTextHelpFormatter
        )
        pdb_subparsers = pdb_parser.add_subparsers(dest="command", required=True, help="PDB处理命令")

        # pdb fetch命令
        fetch_parser = pdb_subparsers.add_parser(
            "fetch", 
            help="PDB文件下载、配体提取和分类管理，以及基于小分子SMILES的结构相似性计算。",
            description="""
            根据提供的蛋白质名称（或基因名称），搜索对应的人类蛋白质UniProt ID，
            然后自动下载相关PDB文件 (mmCIF格式)，提取蛋白质中的配体信息，
            生成报告，并根据是否含有配体将PDB文件分类。
            同时，下载配体化学信息并提取配体坐标。
            使用--smiles参数可根据输入的小分子SMILES计算结构相似性并排序PDB结构。
            """
        )
        fetch_parser.add_argument("protein_name", type=str, help="蛋白质名称或基因名称，例如: BRCA1。")
        fetch_parser.add_argument(
            "--output_dir", 
            type=str, 
            default=None, 
            help="输出目录，用于保存PDB文件和分析结果 (默认: result/pdb_output/蛋白质名_YYYYMMDD_HHMMSS)。")
        fetch_parser.add_argument(
            "--smiles",
            type=str,
            default=None,
            help="输入小分子的SMILES字符串，用于计算与PDB结构中配体的结构相似性并排序PDB结构。")

        # pdb analyze命令
        analyze_parser = pdb_subparsers.add_parser(
            "analyze", 
            help="对指定文件夹下的PDB或CIF文件进行口袋分析，列出配体周围4.5埃内的氨基酸残基。",
            description="""
            对指定文件夹下的PDB或CIF文件中的小分子配体进行口袋分析，
            列出配体距离4.5埃米内的氨基酸残基，
            并将结果写成Markdown文件。
            """
        )
        analyze_parser.add_argument("folder_path", type=str, help="包含PDB或CIF文件的文件夹路径。")
        analyze_parser.add_argument(
            "--output_dir", 
            type=str, 
            default=None, 
            help="输出目录，用于保存分析结果 (默认: 与输入文件夹相同)。")

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
        
        parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.2.0')

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
                if args.command == "fetch":
                    uniprot_id = self._get_human_uniprot_id_and_handle_error(args.protein_name, "PDB处理")
                    if not uniprot_id:
                        return

                    # 设置基于result/的存储路径
                    if args.output_dir is None:
                        output_dir, _ = CommonUtils.get_output_dir("result", "pdb_output", args.protein_name)
                    else:
                        output_dir = args.output_dir

                    self.pdb_processor.process(uniprot_id=uniprot_id, output_dir=output_dir, user_smiles=args.smiles)
                elif args.command == "analyze":
                    # 设置输出目录
                    if args.output_dir is None:
                        output_dir = args.folder_path
                    else:
                        output_dir = args.output_dir

                    self.pdb_processor.analyze_pockets(args.folder_path, output_dir)
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
