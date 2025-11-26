#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
日志管理模块
"""
import sys
from datetime import datetime


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
        print(f"INFO: {info_msg}")
