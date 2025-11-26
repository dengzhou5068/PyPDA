#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
UniProt API交互模块
"""
import json
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from config.config_manager import ConfigManager
from logger.logger import Logger


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
            'User-Agent': 'PyPDA/0.2.0 (https://gitee.com/coding_playground/py-pda; dengzho5068@foxmail.com)'
        })
        return session

    def search_uniprot_by_name(self, name: str, organism_id: int = 9606) -> str:
        """
        根据蛋白质名称或基因名称搜索UniProt，并限定种属。
        Args:
            name: 蛋白质或基因名称。
            organism_id: 限制的物种ID (默认9606为人类)。
        Returns:
            找到的第一个UniProt ID，或None。
        """
        # 定义多种查询策略，从最精确到最宽松
        query_strategies = [
            f"gene_exact:{name} AND organism_id:{organism_id}",  # 精确匹配基因名
            f"gene:{name} AND organism_id:{organism_id}",         # 宽松匹配基因名
            f"name:{name} AND organism_id:{organism_id}",         # 匹配蛋白质名称
            f"{name} AND organism_id:{organism_id}"               # 一般搜索，匹配任何字段
        ]
        
        search_url = f"{self.api_base_url}search"
        
        for i, query in enumerate(query_strategies):
            params = {
                "query": query,
                "format": "json",
                "fields": "accession,protein_name,gene_names",
                "size": 5  # 获取前5个结果进行筛选
            }
            
            try:
                print(f"尝试查询策略 {i+1} 用于 '{name}' (Taxon ID: {organism_id}): {query}")
                response = self.session.get(search_url, params=params, timeout=30)
                response.raise_for_status()
                data = response.json()
                
                if data and 'results' in data and data['results']:
                    # 遍历结果，寻找最匹配的条目
                    for result in data['results']:
                        accession = result.get('primaryAccession')
                        if accession:
                            # 检查基因名或蛋白质名是否匹配
                            gene_names = result.get('genes', [{}])[0].get('geneName', {}).get('value', '').lower()
                            protein_names = result.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', '').lower()
                            name_lower = name.lower()
                            
                            # 如果找到精确匹配或者包含匹配，返回该UniProt ID
                            if name_lower == gene_names or name_lower in gene_names or name_lower in protein_names:
                                print(f"为 '{name}' 找到匹配的 UniProt ID: {accession}")
                                return accession
                    
                    # 如果没有找到精确匹配，但有结果，返回第一个结果
                    accession = data['results'][0].get('primaryAccession')
                    if accession:
                        print(f"为 '{name}' 找到 UniProt ID (使用宽松匹配): {accession}")
                        return accession
            except Exception as e:
                print(f"警告：查询策略 {i+1} 失败: {e}，尝试下一个查询策略")
                continue
        
        # 所有策略都失败后记录错误
        self.logger.log_error(f"无法为 '{name}' (Taxon ID: {organism_id}) 获取UniProt ID，已尝试所有查询策略。", self.config.error_log_path)
        return None

    def get_uniprot_data(self, accession: str) -> dict:
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
            self.logger.log_error(f"网络连接错误 (UniProt API - {url})")
        except requests.exceptions.Timeout:
            self.logger.log_error(f"请求超时 (UniProt API - {url})")
        except json.JSONDecodeError:
            self.logger.log_error(f"UniProt API返回无效JSON (UniProt API - {url}): {response.text[:200]}...")
        except Exception as e:
            self.logger.log_error(f"获取UniProt数据失败: {e} (UniProt API - {url})")
        return None
