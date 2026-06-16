#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
OpenTargets命令处理模块
"""

import csv
import json
from pathlib import Path
from typing import Dict, Any, Optional, List

from config.config_manager import ConfigManager
from logger.logger import Logger
from utils.common_utils import CommonUtils

from .opentargets_api import OpenTargetsAPI


class OpenTargetsProcessor:
    """OpenTargets命令处理类，负责执行各种OpenTargets查询命令"""

    def __init__(self, config: ConfigManager, logger: Logger) -> None:
        """初始化OpenTargetsProcessor实例

        Args:
            config: 配置管理器实例
            logger: 日志记录器实例
        """
        self.config = config
        self.logger = logger
        self.ot_api = OpenTargetsAPI(config, logger)

    def process_disease_drugs(
        self,
        disease_name: Optional[str] = None,
        disease_id: Optional[str] = None,
        limit: int = 100,
        output_dir: Optional[str] = None,
        output_format: str = 'json'
    ) -> None:
        """查询疾病相关药物

        Args:
            disease_name: 疾病名称（与disease_id二选一）
            disease_id: 疾病EFO ID（与disease_name二选一）
            limit: 返回结果数量限制
            output_dir: 输出目录
            output_format: 输出格式（json/csv/all）
        """
        try:
            if disease_id:
                selected_disease_id = disease_id
            elif disease_name:
                selected_disease_id = self.ot_api.get_disease_efo_id(disease_name, interactive=True)
                if not selected_disease_id:
                    return
            else:
                self.logger.log_error("必须提供disease_name或disease_id参数", self.config.error_log_path)
                return

            self.logger.log_info(f"开始查询疾病 {selected_disease_id} 的药物信息", self.config.info_log_path)
            data = self.ot_api.query_disease_drugs(selected_disease_id, limit)
            if not data:
                return

            if output_dir is None:
                output_dir, _ = CommonUtils.get_output_dir("result", "opentargets", selected_disease_id)
            else:
                Path(output_dir).mkdir(parents=True, exist_ok=True)

            json_filename = Path(output_dir) / f"{selected_disease_id}_drugs_{limit}.json"
            with open(json_filename, 'w', encoding='utf-8') as f:
                json.dump(data, f, ensure_ascii=False, indent=4)
            print(f"JSON结果已保存到: {json_filename}")

            if output_format == 'csv' or output_format == 'all':
                csv_filename = json_filename.with_suffix('.csv')
                self._write_disease_drugs_csv(data, csv_filename)
                print(f"CSV结果已保存到: {csv_filename}")

        except Exception as e:
            error_msg = f"处理疾病药物查询失败: {str(e)}"
            self.logger.log_error(error_msg, self.config.error_log_path)

    def process_target_associations(
        self,
        gene_name: str,
        limit: int = 100,
        output_dir: Optional[str] = None,
        output_format: str = 'json'
    ) -> None:
        """查询基因关联疾病

        Args:
            gene_name: 基因名称
            limit: 返回结果数量限制
            output_dir: 输出目录
            output_format: 输出格式（json/csv/all）
        """
        try:
            gene_id = self.ot_api.get_ensembl_id(gene_name)
            if not gene_id:
                return

            self.logger.log_info(f"开始查询基因 {gene_id} 关联的疾病", self.config.info_log_path)
            data = self.ot_api.query_target_associations(gene_id, limit)
            if not data:
                return

            if output_dir is None:
                output_dir, _ = CommonUtils.get_output_dir("result", "opentargets", gene_id)
            else:
                Path(output_dir).mkdir(parents=True, exist_ok=True)

            json_filename = Path(output_dir) / f"{gene_id}_associations_{limit}.json"
            with open(json_filename, 'w', encoding='utf-8') as f:
                json.dump(data, f, ensure_ascii=False, indent=4)
            print(f"JSON结果已保存到: {json_filename}")

            if output_format == 'csv' or output_format == 'all':
                csv_filename = json_filename.with_suffix('.csv')
                self._write_target_associations_csv(data, csv_filename)
                print(f"CSV结果已保存到: {csv_filename}")

        except Exception as e:
            error_msg = f"处理基因关联疾病查询失败: {str(e)}"
            self.logger.log_error(error_msg, self.config.error_log_path)

    def process_target_drugs(
        self,
        gene_name: str,
        output_dir: Optional[str] = None,
        output_format: str = 'json'
    ) -> None:
        """查询靶点相关药物

        Args:
            gene_name: 基因名称
            output_dir: 输出目录
            output_format: 输出格式（json/csv/all）
        """
        try:
            gene_id = self.ot_api.get_ensembl_id(gene_name)
            if not gene_id:
                return

            self.logger.log_info(f"开始查询靶点 {gene_id} 的相关药物", self.config.info_log_path)
            data = self.ot_api.query_target_drugs(gene_id)
            if not data:
                return

            if output_dir is None:
                output_dir, _ = CommonUtils.get_output_dir("result", "opentargets", gene_id)
            else:
                Path(output_dir).mkdir(parents=True, exist_ok=True)

            json_filename = Path(output_dir) / f"{gene_id}_drugs.json"
            with open(json_filename, 'w', encoding='utf-8') as f:
                json.dump(data, f, ensure_ascii=False, indent=4)
            print(f"JSON结果已保存到: {json_filename}")

            if output_format == 'csv' or output_format == 'all':
                csv_filename = json_filename.with_suffix('.csv')
                self._write_target_drugs_csv(data, csv_filename)
                print(f"CSV结果已保存到: {csv_filename}")

        except Exception as e:
            error_msg = f"处理靶点药物查询失败: {str(e)}"
            self.logger.log_error(error_msg, self.config.error_log_path)

    def process_disease_associations(
        self,
        disease_name: Optional[str] = None,
        disease_id: Optional[str] = None,
        limit: int = 100,
        output_dir: Optional[str] = None,
        output_format: str = 'json'
    ) -> None:
        """查询疾病关联靶点

        Args:
            disease_name: 疾病名称（与disease_id二选一）
            disease_id: 疾病EFO ID（与disease_name二选一）
            limit: 返回结果数量限制
            output_dir: 输出目录
            output_format: 输出格式（json/csv/all）
        """
        try:
            if disease_id:
                selected_disease_id = disease_id
            elif disease_name:
                selected_disease_id = self.ot_api.get_disease_efo_id(disease_name, interactive=True)
                if not selected_disease_id:
                    return
            else:
                self.logger.log_error("必须提供disease_name或disease_id参数", self.config.error_log_path)
                return

            self.logger.log_info(f"开始查询疾病 {selected_disease_id} 关联的靶点", self.config.info_log_path)
            data = self.ot_api.query_disease_associations(selected_disease_id, limit)
            if not data:
                return

            if output_dir is None:
                output_dir, _ = CommonUtils.get_output_dir("result", "opentargets", selected_disease_id)
            else:
                Path(output_dir).mkdir(parents=True, exist_ok=True)

            json_filename = Path(output_dir) / f"{selected_disease_id}_targets_{limit}.json"
            with open(json_filename, 'w', encoding='utf-8') as f:
                json.dump(data, f, ensure_ascii=False, indent=4)
            print(f"JSON结果已保存到: {json_filename}")

            if output_format == 'csv' or output_format == 'all':
                csv_filename = json_filename.with_suffix('.csv')
                self._write_disease_associations_csv(data, csv_filename)
                print(f"CSV结果已保存到: {csv_filename}")

        except Exception as e:
            error_msg = f"处理疾病关联靶点查询失败: {str(e)}"
            self.logger.log_error(error_msg, self.config.error_log_path)

    def _write_disease_drugs_csv(self, data: Dict[str, Any], csv_path: Path) -> None:
        """写入疾病药物信息的CSV格式"""
        try:
            with open(csv_path, 'w', newline='', encoding='utf-8') as csvfile:
                headers = ["drug.id", "drug.name", "maxClinicalStage", "drugType", "tradeNames"]
                writer = csv.writer(csvfile)
                writer.writerow(headers)

                drugs = []
                if isinstance(data, dict) and 'data' in data and isinstance(data['data'], dict):
                    disease = data['data'].get('disease', {})
                    if isinstance(disease, dict):
                        drug_candidates = disease.get('drugAndClinicalCandidates', {})
                        if isinstance(drug_candidates, dict):
                            drugs = drug_candidates.get('rows', [])

                if not isinstance(drugs, list):
                    drugs = []

                for drug_info in drugs:
                    if not isinstance(drug_info, dict):
                        continue
                    drug = drug_info.get('drug', {})
                    if not isinstance(drug, dict):
                        drug = {}

                    trade_names = drug.get('tradeNames', [])
                    trade_names_str = ', '.join(trade_names) if isinstance(trade_names, list) else ''

                    writer.writerow([
                        drug.get('id', ''),
                        drug.get('name', ''),
                        drug_info.get('maxClinicalStage', ''),
                        drug.get('drugType', ''),
                        trade_names_str
                    ])
        except Exception as e:
            error_msg = f"写入疾病药物CSV失败: {str(e)}"
            self.logger.log_error(error_msg, self.config.error_log_path)

    def _write_target_drugs_csv(self, data: Dict[str, Any], csv_path: Path) -> None:
        """写入靶点药物信息的CSV格式"""
        try:
            with open(csv_path, 'w', newline='', encoding='utf-8') as csvfile:
                headers = ["drug.id", "drug.name", "maxClinicalStage", "drugType", "tradeNames", "diseases"]
                writer = csv.writer(csvfile)
                writer.writerow(headers)

                drugs = []
                if isinstance(data, dict) and 'data' in data and isinstance(data['data'], dict):
                    target = data['data'].get('target', {})
                    if isinstance(target, dict):
                        drug_candidates = target.get('drugAndClinicalCandidates', {})
                        if isinstance(drug_candidates, dict):
                            drugs = drug_candidates.get('rows', [])

                if not isinstance(drugs, list):
                    drugs = []

                for drug_info in drugs:
                    if not isinstance(drug_info, dict):
                        continue
                    drug = drug_info.get('drug', {})
                    if not isinstance(drug, dict):
                        drug = {}

                    trade_names = drug.get('tradeNames', [])
                    trade_names_str = ', '.join(trade_names) if isinstance(trade_names, list) else ''

                    diseases = drug_info.get('diseases', [])
                    if isinstance(diseases, list):
                        disease_names = [d.get('diseaseFromSource', '') for d in diseases if isinstance(d, dict)]
                        diseases_str = ', '.join(disease_names)
                    else:
                        diseases_str = ''

                    writer.writerow([
                        drug.get('id', ''),
                        drug.get('name', ''),
                        drug_info.get('maxClinicalStage', ''),
                        drug.get('drugType', ''),
                        trade_names_str,
                        diseases_str
                    ])
        except Exception as e:
            error_msg = f"写入靶点药物CSV失败: {str(e)}"
            self.logger.log_error(error_msg, self.config.error_log_path)

    def _write_target_associations_csv(self, data: Dict[str, Any], csv_path: Path) -> None:
        """写入目标关联疾病CSV数据"""
        try:
            rows = []
            if isinstance(data, dict) and 'data' in data and isinstance(data['data'], dict):
                target = data['data'].get('target', {})
                if isinstance(target, dict):
                    associated_diseases = target.get('associatedDiseases', {})
                    if isinstance(associated_diseases, dict):
                        rows = associated_diseases.get('rows', [])

            if not isinstance(rows, list):
                rows = []

            all_datasources = set()
            for row in rows:
                if isinstance(row, dict):
                    scores = row.get('datasourceScores', [])
                    if isinstance(scores, list):
                        for score in scores:
                            if isinstance(score, dict):
                                all_datasources.add(score.get('componentId', ''))

            datasource_columns = sorted(all_datasources)
            headers = ["disease.id", "disease.name", "score"] + datasource_columns

            with open(csv_path, 'w', newline='', encoding='utf-8') as csvfile:
                writer = csv.writer(csvfile)
                writer.writerow(headers)

                for row in rows:
                    if not isinstance(row, dict):
                        continue

                    disease = row.get('disease', {})
                    if not isinstance(disease, dict):
                        disease = {}

                    scores = row.get('datasourceScores', [])
                    if not isinstance(scores, list):
                        scores = []

                    score_map = {}
                    for score in scores:
                        if isinstance(score, dict):
                            score_map[score.get('componentId', '')] = str(score.get('score', ''))

                    row_data = [
                        disease.get('id', ''),
                        disease.get('name', ''),
                        str(row.get('score', ''))
                    ] + [score_map.get(ds, '') for ds in datasource_columns]

                    writer.writerow(row_data)
        except Exception as e:
            error_msg = f"写入目标关联疾病CSV失败: {str(e)}"
            self.logger.log_error(error_msg, self.config.error_log_path)

    def _write_disease_associations_csv(self, data: Dict[str, Any], csv_path: Path) -> None:
        """写入疾病关联靶点CSV数据"""
        try:
            rows = []
            if isinstance(data, dict) and 'data' in data and isinstance(data['data'], dict):
                disease = data['data'].get('disease', {})
                if isinstance(disease, dict):
                    associated_targets = disease.get('associatedTargets', {})
                    if isinstance(associated_targets, dict):
                        rows = associated_targets.get('rows', [])

            if not isinstance(rows, list):
                rows = []

            all_datasources = set()
            for row in rows:
                if isinstance(row, dict):
                    scores = row.get('datasourceScores', [])
                    if isinstance(scores, list):
                        for score in scores:
                            if isinstance(score, dict):
                                all_datasources.add(score.get('componentId', ''))

            datasource_columns = sorted(all_datasources)
            headers = ["target.id", "target.approvedSymbol", "target.approvedName", "score"] + datasource_columns

            with open(csv_path, 'w', newline='', encoding='utf-8') as csvfile:
                writer = csv.writer(csvfile)
                writer.writerow(headers)

                for row in rows:
                    if not isinstance(row, dict):
                        continue

                    target = row.get('target', {})
                    if not isinstance(target, dict):
                        target = {}

                    scores = row.get('datasourceScores', [])
                    if not isinstance(scores, list):
                        scores = []

                    score_map = {}
                    for score in scores:
                        if isinstance(score, dict):
                            score_map[score.get('componentId', '')] = str(score.get('score', ''))

                    row_data = [
                        target.get('id', ''),
                        target.get('approvedSymbol', ''),
                        target.get('approvedName', ''),
                        str(row.get('score', ''))
                    ] + [score_map.get(ds, '') for ds in datasource_columns]

                    writer.writerow(row_data)
        except Exception as e:
            error_msg = f"写入疾病关联靶点CSV失败: {str(e)}"
            self.logger.log_error(error_msg, self.config.error_log_path)
