#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
配置管理模块
"""
import os
from typing import Dict


class ConfigManager:
    """配置文件管理类，负责读取和管理配置信息"""
    def __init__(self) -> None:
        """初始化ConfigManager实例
        
        设置默认配置项，包括命名空间、日志文件路径和API基础URL
        """
        self.ns: Dict[str, str] = {'uniprot': 'http://uniprot.org/uniprot'}
        self.error_log_path: str = os.path.join(os.getcwd(), 'error.txt')
        self.info_log_path: str = os.path.join(os.getcwd(), 'info.txt')
        self.uniprot_api_base_url: str = 'https://rest.uniprot.org/uniprotkb/'
        self.opentargets_api_base_url: str = 'https://api.platform.opentargets.org/api/v4/graphql'
        self.ensembl_api_url: str = 'https://rest.ensembl.org/lookup/symbol/homo_sapiens/'
        self.ols_url: str = 'https://www.ebi.ac.uk/ols/api/search'
        self.kegg_api_base_url: str = 'https://rest.kegg.jp/'
