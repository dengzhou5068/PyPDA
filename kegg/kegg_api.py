#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
KEGG API交互模块
"""
import json
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from typing import Dict, List, Optional, Any, Tuple

from config.config_manager import ConfigManager
from logger.logger import Logger


class KEGGAPI:
    """KEGG API交互类，负责从KEGG API获取数据"""

    def __init__(self, config: ConfigManager, logger: Logger) -> None:
        """初始化KEGGAPI实例

        Args:
            config: 配置管理器实例
            logger: 日志记录器实例
        """
        self.config = config
        self.logger = logger
        self.api_base_url = config.kegg_api_base_url
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
            'User-Agent': (
                'PyPDA/0.6.0 '
                '(https://github.com/dengzhou5068/PyPDA; '
                'dengzho5068@foxmail.com)'
            )
        })
        return session

    def search_gene_by_name(self, gene_name: str, organism: str = 'hsa') -> Optional[str]:
        """
        根据基因名称搜索KEGG基因ID

        Args:
            gene_name: 基因名称
            organism: 生物体代码，默认为hsa（人类）

        Returns:
            KEGG基因ID（如hsa:10000），或None
        """
        search_urls = [
            f"{self.api_base_url}find/genes/{organism}:{gene_name}",
            f"{self.api_base_url}find/genes/{gene_name}",
            f"{self.api_base_url}list/genes/{organism}"
        ]

        for url in search_urls:
            try:
                response = self.session.get(url, timeout=30)
                response.raise_for_status()
                lines = response.text.strip().split('\n')
                if lines:
                    for line in lines:
                        if line:
                            parts = line.split('\t')
                            if len(parts) >= 2:
                                gene_id = parts[0]
                                gene_symbol = parts[1]
                                if gene_symbol.lower() == gene_name.lower() or gene_name.lower() in gene_symbol.lower():
                                    print(f"为基因 '{gene_name}' 找到KEGG基因ID: {gene_id}")
                                    return gene_id
                    # 如果没有找到精确匹配，返回第一个结果
                    first_line = lines[0]
                    gene_id = first_line.split('\t')[0]
                    print(f"为基因 '{gene_name}' 找到KEGG基因ID (宽松匹配): {gene_id}")
                    return gene_id
            except Exception:
                continue

        error_msg = f"无法为基因 '{gene_name}' 找到KEGG基因ID"
        self.logger.log_error(error_msg, self.config.error_log_path)
        return None

    def convert_uniprot_to_kegg(self, uniprot_id: str, organism: str = 'hsa') -> Optional[str]:
        """
        使用KEGG conv操作将UniProt ID转换为KEGG基因ID

        Args:
            uniprot_id: UniProt蛋白质ID（如P04637）
            organism: 生物体代码，默认为hsa（人类）

        Returns:
            KEGG基因ID（如hsa:7159），或None
        """
        url = f"{self.api_base_url}conv/{organism}/uniprot:{uniprot_id}"
        try:
            response = self.session.get(url, timeout=30)
            response.raise_for_status()
            lines = response.text.strip().split('\n')
            if lines and lines[0]:
                parts = lines[0].split('\t')
                if len(parts) >= 2:
                    kegg_id = parts[1]
                    print(f"将UniProt ID '{uniprot_id}' 转换为KEGG基因ID: {kegg_id}")
                    return kegg_id
            return None
        except requests.exceptions.HTTPError as e:
            error_msg = f"HTTP请求错误 (KEGG conv - {url}): "
            error_msg += f"{e.response.status_code} - {e.response.text}"
            self.logger.log_error(error_msg, self.config.error_log_path)
        except Exception as e:
            error_msg = f"UniProt ID转换失败: {e} (KEGG conv - {url})"
            self.logger.log_error(error_msg, self.config.error_log_path)
        return None

    def get_pathways_by_gene(self, gene_id: str) -> List[Tuple[str, str]]:
        """
        获取基因参与的信号通路列表

        Args:
            gene_id: KEGG基因ID（如hsa:10000）

        Returns:
            信号通路ID和名称的元组列表
        """
        pathways = []
        url = f"{self.api_base_url}link/pathway/{gene_id}"
        try:
            response = self.session.get(url, timeout=30)
            response.raise_for_status()
            lines = response.text.strip().split('\n')
            for line in lines:
                if line:
                    parts = line.split('\t')
                    if len(parts) >= 2:
                        pathway_id = parts[1]
                        pathway_name = self._get_pathway_name(pathway_id)
                        pathways.append((pathway_id, pathway_name))
            print(f"基因 {gene_id} 参与 {len(pathways)} 个信号通路")
            return pathways
        except requests.exceptions.HTTPError as e:
            error_msg = f"HTTP请求错误 (KEGG API - {url}): "
            error_msg += f"{e.response.status_code} - {e.response.text}"
            self.logger.log_error(error_msg, self.config.error_log_path)
        except Exception as e:
            error_msg = f"获取信号通路失败: {e} (KEGG API - {url})"
            self.logger.log_error(error_msg, self.config.error_log_path)
        return pathways

    def _get_pathway_name(self, pathway_id: str) -> str:
        """
        获取信号通路的名称

        Args:
            pathway_id: 信号通路ID（如hsa04110）

        Returns:
            信号通路名称
        """
        url = f"{self.api_base_url}get/{pathway_id}"
        try:
            response = self.session.get(url, timeout=30)
            response.raise_for_status()
            lines = response.text.strip().split('\n')
            for line in lines:
                if line.startswith('NAME'):
                    return line.split(' ', 1)[1].strip()
            return pathway_id
        except Exception:
            return pathway_id

    def download_pathway_image(self, pathway_id: str, output_file: str) -> bool:
        """
        下载KEGG信号通路图

        Args:
            pathway_id: 信号通路ID（如path:hsa04110或hsa04110）
            output_file: 输出文件路径（应包含.png扩展名）

        Returns:
            如果下载成功返回True，否则返回False
        """
        clean_pathway_id = pathway_id.replace('path:', '')
        urls = [
            f"https://rest.kegg.jp/get/{clean_pathway_id}/image",
            f"https://www.kegg.jp/kegg/pathway/hsa/{clean_pathway_id}.png",
            f"https://www.genome.jp/kegg/pathway/{clean_pathway_id}/{clean_pathway_id}.png"
        ]

        for url in urls:
            try:
                response = self.session.get(url, timeout=60)
                response.raise_for_status()
                with open(output_file, 'wb') as f:
                    f.write(response.content)
                print(f"信号通路图已保存至 {output_file}")
                return True
            except requests.exceptions.HTTPError:
                continue
            except Exception as e:
                continue

        error_msg = f"下载信号通路图失败: 尝试了所有URL都无法获取 {clean_pathway_id}"
        self.logger.log_error(error_msg, self.config.error_log_path)
        return False

    def get_gene_info(self, gene_id: str) -> Optional[Dict[str, Any]]:
        """
        获取基因详细信息

        Args:
            gene_id: KEGG基因ID（如hsa:10000）

        Returns:
            基因信息字典，或None
        """
        url = f"{self.api_base_url}get/{gene_id}"
        try:
            response = self.session.get(url, timeout=30)
            response.raise_for_status()
            return self._parse_kegg_entry(response.text)
        except requests.exceptions.HTTPError as e:
            error_msg = f"HTTP请求错误 (KEGG API - {url}): "
            error_msg += f"{e.response.status_code} - {e.response.text}"
            self.logger.log_error(error_msg, self.config.error_log_path)
        except Exception as e:
            error_msg = f"获取基因信息失败: {e} (KEGG API - {url})"
            self.logger.log_error(error_msg, self.config.error_log_path)
        return None

    def _parse_kegg_entry(self, text: str) -> Dict[str, Any]:
        """
        解析KEGG条目文本

        Args:
            text: KEGG条目文本

        Returns:
            解析后的字典
        """
        result = {}
        current_section = None
        current_value = []

        for line in text.split('\n'):
            if not line:
                continue

            if line.startswith('///'):
                break

            if line[0].isspace():
                if current_section:
                    current_value.append(line.strip())
            else:
                if current_section:
                    result[current_section] = '\n'.join(current_value)
                parts = line.split(' ', 1)
                current_section = parts[0]
                current_value = [parts[1].strip()] if len(parts) > 1 else []

        if current_section:
            result[current_section] = '\n'.join(current_value)

        return result

    def search_pathway_by_name(self, pathway_name: str) -> List[Tuple[str, str]]:
        """
        根据名称搜索信号通路

        Args:
            pathway_name: 信号通路名称

        Returns:
            信号通路ID和名称的元组列表
        """
        pathways = []
        url = f"{self.api_base_url}find/pathway/{pathway_name}"
        try:
            response = self.session.get(url, timeout=30)
            response.raise_for_status()
            lines = response.text.strip().split('\n')
            for line in lines:
                if line:
                    parts = line.split('\t')
                    if len(parts) >= 2:
                        pathways.append((parts[0], parts[1]))
            return pathways
        except requests.exceptions.HTTPError as e:
            error_msg = f"HTTP请求错误 (KEGG API - {url}): "
            error_msg += f"{e.response.status_code} - {e.response.text}"
            self.logger.log_error(error_msg, self.config.error_log_path)
        except Exception as e:
            error_msg = f"搜索信号通路失败: {e} (KEGG API - {url})"
            self.logger.log_error(error_msg, self.config.error_log_path)
        return pathways
