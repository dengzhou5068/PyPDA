#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
配置管理模块
"""
import os


class ConfigManager:
    """配置文件管理类，负责读取和管理配置信息"""
    def __init__(self):
        self.ns = {'uniprot': 'http://uniprot.org/uniprot'}
        self.error_log_path = os.path.join(os.getcwd(), 'error.txt')
        self.info_log_path = os.path.join(os.getcwd(), 'info.txt')
        self.uniprot_api_base_url = 'https://rest.uniprot.org/uniprotkb/'
