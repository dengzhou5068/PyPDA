#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
通用工具模块
"""
import json
import os
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
from pathlib import Path
from typing import List, Dict, Tuple, Any, Union

from Bio import SeqIO
from tqdm import tqdm

from logger.logger import Logger


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
            Logger.log_error(f"FASTA文件未找到: {file_path}", "error.txt")
            raise
        except Exception as e:
            Logger.log_error(f"读取FASTA文件失败: {file_path} - {str(e)}", "error.txt")
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
            Logger.log_error(f"保存序列到文件失败: {output_file} - {str(e)}", "error.txt")
            raise

    @staticmethod
    def parallel_executor(func: callable, items: List[Any], max_workers: Union[int, None] = None, description: str = "Processing") -> None:
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
            Logger.log_error(f"保存JSON文件失败: {json_filename} - {str(e)}", "error.txt")
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
            Logger.log_error(f"文件 {file_path} 不存在", "error.txt")
            raise FileNotFoundError(f"文件 {file_path} 不存在")

        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                return json.load(f)
        except json.JSONDecodeError as e:
            Logger.log_error(f"JSON文件解析错误: {file_path} - {str(e)}", "error.txt")
            raise
        except IOError as e:
            Logger.log_error(f"读取文件失败: {file_path} - {str(e)}", "error.txt")
            raise

    @staticmethod
    def get_output_dir(base_dir: str, tool_name: str, input_name: str = "output", timestamp: str = None) -> Tuple[str, str]:
        """根据工具名称和输入名称生成带时间戳的输出目录

        Args:
            base_dir: 基础目录，如 "result"
            tool_name: 工具名称，如 "protein_sequences", "pdb_output", "uniprot_reports"
            input_name: 输入名称，如基因名或蛋白质名
            timestamp: 时间戳，如果为None则使用当前时间

        Returns:
            完整的输出目录路径和时间戳
        """
        if timestamp is None:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_dir = Path(base_dir) / tool_name / f"{input_name}_{timestamp}"
        output_dir.mkdir(parents=True, exist_ok=True)
        return str(output_dir), timestamp
