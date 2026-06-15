#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
OpenTargets API交互模块
"""
import json
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from typing import Dict, Any, Optional, List

from config.config_manager import ConfigManager
from logger.logger import Logger


class OpenTargetsAPI:
    """OpenTargets API交互类，负责从OpenTargets API获取数据"""

    def __init__(self, config: ConfigManager, logger: Logger) -> None:
        """初始化OpenTargetsAPI实例

        Args:
            config: 配置管理器实例
            logger: 日志记录器实例
        """
        self.config = config
        self.logger = logger
        self.opentargets_api_base_url = config.opentargets_api_base_url
        self.ensembl_api_url = config.ensembl_api_url
        self.ols_url = config.ols_url
        self.session = self._setup_session()

    def _setup_session(self) -> requests.Session:
        """设置带有重试策略的requests会话"""
        retry_strategy = Retry(
            total=5,
            backoff_factor=1,
            status_forcelist=[429, 500, 502, 503, 504],
            allowed_methods=["HEAD", "GET", "POST", "OPTIONS"]
        )
        adapter = HTTPAdapter(max_retries=retry_strategy)
        session = requests.Session()
        session.mount('http://', adapter)
        session.mount('https://', adapter)
        session.headers.update({
            'Content-Type': 'application/json',
            'Accept': 'application/json',
            'User-Agent': (
                'PyPDA/0.5.0 '
                '(https://gitee.com/coding_playground/py-pda; '
                'dengzho5068@foxmail.com)'
            )
        })
        return session

    def get_ensembl_id(self, gene_name: str) -> Optional[str]:
        """根据基因名称获取Ensembl ID

        Args:
            gene_name: 基因名称

        Returns:
            Ensembl ID，如果未找到则返回None
        """
        url = f"{self.ensembl_api_url}{gene_name}"
        params = {'content-type': 'application/json'}

        try:
            response = self.session.get(url, params=params, timeout=30)
            response.raise_for_status()
            data = response.json()
            # Ensembl API可能返回列表或字典
            if isinstance(data, list) and len(data) > 0:
                target_id = data[0].get("id")
            elif isinstance(data, dict):
                target_id = data.get("id")
            else:
                target_id = None
                
            if not target_id:
                error_msg = f"未找到基因 '{gene_name}' 的Ensembl ID"
                self.logger.log_error(error_msg, self.config.error_log_path)
                return None
            return target_id
        except requests.exceptions.HTTPError as e:
            error_msg = f"HTTP请求错误 (Ensembl API - {url}): "
            error_msg += f"{e.response.status_code} - {e.response.text}"
            self.logger.log_error(error_msg, self.config.error_log_path)
        except requests.exceptions.ConnectionError:
            self.logger.log_error(f"网络连接错误 (Ensembl API - {url})", self.config.error_log_path)
        except requests.exceptions.Timeout:
            self.logger.log_error(f"请求超时 (Ensembl API - {url})", self.config.error_log_path)
        except json.JSONDecodeError:
            error_msg = f"Ensembl API返回无效JSON (Ensembl API - {url}): "
            error_msg += f"{response.text[:200]}..."
            self.logger.log_error(error_msg, self.config.error_log_path)
        except Exception as e:
            error_msg = f"获取Ensembl ID失败: {e} (Ensembl API - {url})"
            self.logger.log_error(error_msg, self.config.error_log_path)
        return None

    def get_disease_efo_id(self, disease_name: str, interactive: bool = True) -> Optional[str]:
        """根据疾病名称获取EFO ID

        Args:
            disease_name: 疾病名称
            interactive: 是否使用交互式选择（当有多个匹配结果时）

        Returns:
            疾病的EFO ID，如果未找到则返回None
        """
        params = {
            "q": disease_name,
            "ontology": "efo",
            "exact": "true"
        }

        try:
            response = self.session.get(self.ols_url, params=params, timeout=30)
            response.raise_for_status()
            data = response.json()

            if data["response"]["numFound"] == 0:
                error_msg = f"未找到疾病 '{disease_name}' 的EFO ID"
                self.logger.log_error(error_msg, self.config.error_log_path)
                return None

            docs = data["response"]["docs"]
            if not docs or len(docs) == 0:
                error_msg = f"OLS API返回空结果列表 (疾病: {disease_name})"
                self.logger.log_error(error_msg, self.config.error_log_path)
                return None
                
            sorted_disease_hits = sorted(
                docs,
                key=lambda x: x["label"].lower().startswith(disease_name.lower()),
                reverse=True
            )

            if len(sorted_disease_hits) == 1:
                return sorted_disease_hits[0]["iri"].split("/")[-1]
            else:
                if interactive:
                    return self._select_disease_interactive(sorted_disease_hits, disease_name)
                else:
                    self.logger.log_info(
                        f"找到多个匹配的疾病: {disease_name}，使用第一个匹配",
                        self.config.info_log_path
                    )
                    return sorted_disease_hits[0]["iri"].split("/")[-1]
        except requests.exceptions.HTTPError as e:
            error_msg = f"HTTP请求错误 (OLS API - {self.ols_url}): "
            error_msg += f"{e.response.status_code} - {e.response.text}"
            self.logger.log_error(error_msg, self.config.error_log_path)
        except requests.exceptions.ConnectionError:
            self.logger.log_error(f"网络连接错误 (OLS API - {self.ols_url})", self.config.error_log_path)
        except requests.exceptions.Timeout:
            self.logger.log_error(f"请求超时 (OLS API - {self.ols_url})", self.config.error_log_path)
        except json.JSONDecodeError:
            error_msg = f"OLS API返回无效JSON (OLS API - {self.ols_url}): "
            error_msg += f"{response.text[:200]}..."
            self.logger.log_error(error_msg, self.config.error_log_path)
        except Exception as e:
            error_msg = f"获取EFO ID失败: {e} (OLS API - {self.ols_url})"
            self.logger.log_error(error_msg, self.config.error_log_path)
        return None

    def _select_disease_interactive(self, docs: List[Dict[str, Any]], disease_name: str) -> Optional[str]:
        """交互式选择疾病"""
        print(f"\n找到 {len(docs)} 个匹配的疾病，请选择要查询的疾病：")
        print("-" * 60)

        for index, doc in enumerate(docs, start=1):
            efo_id = doc["iri"].split("/")[-1]
            label = doc.get('label', '未知')
            description = doc.get('description', ['无描述'])[0] if isinstance(doc.get('description'), list) else doc.get('description', '无描述')
            print(f"{index}. {label}")
            print(f"   EFO ID: {efo_id}")
            if description and description != '无描述':
                print(f"   描述: {description[:100]}{'...' if len(description) > 100 else ''}")
            print()

        while True:
            try:
                choice = input(
                    f"请输入要选择的序号 (1-{len(docs)})，或按 q 退出: "
                ).strip()
                if choice.lower() == 'q':
                    self.logger.log_info("用户取消操作", self.config.info_log_path)
                    return None

                choice_num = int(choice)
                if 1 <= choice_num <= len(docs):
                    selected_doc = docs[choice_num - 1]
                    selected_id = selected_doc["iri"].split("/")[-1]
                    selected_label = selected_doc.get('label', '未知')
                    print(f"已选择: {selected_label} ({selected_id})")
                    return selected_id
                else:
                    print(f"请输入 1-{len(docs)} 之间的数字")
            except ValueError:
                print("请输入有效的数字或 q 退出")

    def send_graphql_request(self, query: str, variables: Dict[str, Any]) -> Optional[Dict[str, Any]]:
        """发送GraphQL请求

        Args:
            query: GraphQL查询字符串
            variables: 查询变量

        Returns:
            API返回的JSON数据，如果失败则返回None
        """
        payload = {"query": query, "variables": variables}

        try:
            response = self.session.post(self.opentargets_api_base_url, json=payload, timeout=60)
            response.raise_for_status()
            return response.json()
        except requests.exceptions.HTTPError as e:
            error_msg = f"HTTP请求错误 (OpenTargets API - {self.opentargets_api_base_url}): "
            error_msg += f"{e.response.status_code} - {e.response.text}"
            self.logger.log_error(error_msg, self.config.error_log_path)
        except requests.exceptions.ConnectionError:
            self.logger.log_error(f"网络连接错误 (OpenTargets API - {self.opentargets_api_base_url})", self.config.error_log_path)
        except requests.exceptions.Timeout:
            self.logger.log_error(f"请求超时 (OpenTargets API - {self.opentargets_api_base_url})", self.config.error_log_path)
        except json.JSONDecodeError:
            error_msg = f"OpenTargets API返回无效JSON: {response.text[:200]}..."
            self.logger.log_error(error_msg, self.config.error_log_path)
        except Exception as e:
            error_msg = f"发送GraphQL请求失败: {e}"
            self.logger.log_error(error_msg, self.config.error_log_path)
        return None

    def query_disease_drugs(self, disease_id: str, limit: int = 100) -> Optional[Dict[str, Any]]:
        """查询特定疾病的药物信息

        Args:
            disease_id: 疾病的EFO ID
            limit: 返回结果数量限制

        Returns:
            药物信息数据，如果失败则返回None
        """
        query = """
        query DiseaseDrugs($diseaseId: String!, $size: Int!) {
          disease(efoId: $diseaseId) {
            id
            name
            associatedTargets(
              page: {index: 0, size: $size}
              orderByScore: "score"
              enableIndirect: true
            ) {
              count
              rows {
                target {
                  id
                  approvedSymbol
                  approvedName
                  tractability {
                    modality
                    label
                  }
                }
                score
              }
            }
          }
        }
        """

        variables = {
            "diseaseId": disease_id,
            "size": limit
        }
        return self.send_graphql_request(query, variables)

    def query_target_associations(self, gene_id: str, limit: int = 100) -> Optional[Dict[str, Any]]:
        """查询特定基因关联的疾病

        Args:
            gene_id: 基因的Ensembl ID
            limit: 返回结果数量限制

        Returns:
            疾病关联数据，如果失败则返回None
        """
        query = """
        query TargetAssociations($id: String!, $size: Int!) {
          target(ensemblId: $id) {
            id
            approvedSymbol
            associatedDiseases(
              page: { index: 0, size: $size }
              orderByScore: "score"
              enableIndirect: false
            ) {
              count
              rows {
                disease { id name }
                score
                datasourceScores { componentId: id score }
              }
            }
          }
        }
        """

        variables = {
            "id": gene_id,
            "size": limit
        }
        return self.send_graphql_request(query, variables)

    def query_disease_associations(self, disease_id: str, limit: int = 100) -> Optional[Dict[str, Any]]:
        """查询特定疾病关联的靶点

        Args:
            disease_id: 疾病的EFO ID
            limit: 返回结果数量限制

        Returns:
            靶点关联数据，如果失败则返回None
        """
        query = """
        query DiseaseAssociations($id: String!, $size: Int!) {
          disease(efoId: $id) {
            id
            name
            associatedTargets(
              page: { index: 0, size: $size }
              orderByScore: "score"
              enableIndirect: true
            ) {
              count
              rows {
                target {
                  id
                  approvedSymbol
                  approvedName
                }
                score
                datasourceScores { componentId: id score }
              }
            }
          }
        }
        """

        variables = {
            "id": disease_id,
            "size": limit
        }
        return self.send_graphql_request(query, variables)

    def query_target_drugs(self, gene_id: str, limit: int = 100) -> Optional[Dict[str, Any]]:
        """查询特定基因关联的药物信息

        Args:
            gene_id: 基因的Ensembl ID
            limit: 返回结果数量限制

        Returns:
            药物信息数据，如果失败则返回None
        """
        query = """
        query TargetDrugs($id: String!) {
          target(ensemblId: $id) {
            id
            approvedSymbol
            approvedName
            tractability {
              modality
              label
            }
            associatedDiseases(page: {index: 0, size: 5}) {
              rows {
                disease {
                  id
                  name
                }
              }
            }
          }
        }
        """

        variables = {
            "id": gene_id
        }
        return self.send_graphql_request(query, variables)
