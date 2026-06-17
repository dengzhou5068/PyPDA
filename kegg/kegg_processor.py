#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
KEGG处理器模块
"""
import argparse
import csv
import json
from pathlib import Path
from typing import List, Tuple, Optional

from config.config_manager import ConfigManager
from logger.logger import Logger
from kegg.kegg_api import KEGGAPI
from utils.common_utils import CommonUtils


class KEGGProcessor:
    """KEGG处理器类，负责处理KEGG相关命令"""

    def __init__(self, config: ConfigManager, logger: Logger, uniprot_api=None):
        """初始化KEGGProcessor实例

        Args:
            config: 配置管理器实例
            logger: 日志记录器实例
            uniprot_api: UniProtAPI实例，用于将基因名称转换为UniProt ID
        """
        self.config = config
        self.logger = logger
        self.kegg_api = KEGGAPI(config, logger)
        self.uniprot_api = uniprot_api

    def _save_json_with_csv(self, data: List[dict], json_file: Path) -> None:
        """
        同时保存JSON和CSV文件

        Args:
            data: 要保存的数据列表
            json_file: JSON输出文件路径
        """
        with open(json_file, 'w', encoding='utf-8') as f:
            json.dump(data, f, ensure_ascii=False, indent=4)
        print(f"JSON文件已保存至 {json_file}")

        csv_file = json_file.with_suffix('.csv')
        if data:
            headers = data[0].keys()
            with open(csv_file, 'w', encoding='utf-8-sig', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(headers)
                for row in data:
                    writer.writerow(row.values())
            print(f"CSV文件已保存至 {csv_file}")

    def process_gene_pathways(self, gene_name: str, output_dir: Optional[str] = None) -> None:
        """
        根据基因名称查询参与的信号通路并下载通路图
        流程：基因名称 -> UniProt ID -> KEGG基因ID -> 信号通路

        Args:
            gene_name: 基因名称
            output_dir: 输出目录
        """
        if output_dir is None:
            output_dir, _ = CommonUtils.get_output_dir("result", "kegg_output", gene_name)
        else:
            output_dir = str(CommonUtils.ensure_output_dir(output_dir))

        kegg_gene_id = None

        if self.uniprot_api:
            print(f"步骤1: 根据基因名称 '{gene_name}' 查询UniProt ID")
            uniprot_id = self.uniprot_api.search_uniprot_by_name(gene_name)
            if uniprot_id:
                print(f"步骤2: 将UniProt ID '{uniprot_id}' 转换为KEGG基因ID")
                kegg_gene_id = self.kegg_api.convert_uniprot_to_kegg(uniprot_id)

        if not kegg_gene_id:
            print(f"步骤1/2失败，尝试直接搜索KEGG基因ID")
            kegg_gene_id = self.kegg_api.search_gene_by_name(gene_name)

        if not kegg_gene_id:
            print(f"无法为基因 '{gene_name}' 找到KEGG基因ID")
            return

        pathways = self.kegg_api.get_pathways_by_gene(kegg_gene_id)
        if not pathways:
            print(f"基因 '{gene_name}' (KEGG ID: {kegg_gene_id}) 未找到参与的信号通路")
            return

        print(f"基因 '{gene_name}' (KEGG ID: {kegg_gene_id}) 参与以下信号通路:")
        for pathway_id, pathway_name in pathways:
            print(f"  - {pathway_id}: {pathway_name}")

        pathways_data = []
        for pathway_id, pathway_name in pathways:
            pathways_data.append({
                "pathway_id": pathway_id,
                "pathway_name": pathway_name
            })

        json_output_file = Path(output_dir) / f"{gene_name}_pathways.json"
        self._save_json_with_csv(pathways_data, json_output_file)

        pathways_dir = Path(output_dir) / "pathway_images"
        pathways_dir.mkdir(parents=True, exist_ok=True)

        for pathway_id, pathway_name in pathways:
            safe_pathway_id = pathway_id.replace(':', '_')
            image_file = pathways_dir / f"{safe_pathway_id}.png"
            success = self.kegg_api.download_pathway_image(pathway_id, str(image_file))
            if success:
                print(f"成功下载信号通路图: {pathway_name} -> {image_file}")
            else:
                print(f"下载信号通路图失败: {pathway_name}")

    def process_download_pathway(self, pathway_id: str, output_dir: Optional[str] = None) -> None:
        """
        根据信号通路ID下载信号通路图

        Args:
            pathway_id: 信号通路ID（如hsa04110）
            output_dir: 输出目录
        """
        if output_dir is None:
            output_dir, _ = CommonUtils.get_output_dir("result", "kegg_output", pathway_id)
        else:
            output_dir = str(CommonUtils.ensure_output_dir(output_dir))

        image_file = Path(output_dir) / f"{pathway_id}.png"
        success = self.kegg_api.download_pathway_image(pathway_id, str(image_file))
        if success:
            print(f"成功下载信号通路图至 {image_file}")
        else:
            print(f"下载信号通路图失败: {pathway_id}")

    def process_search_pathway(self, pathway_name: str, output_dir: Optional[str] = None) -> None:
        """
        根据名称搜索信号通路

        Args:
            pathway_name: 信号通路名称
            output_dir: 输出目录
        """
        pathways = self.kegg_api.search_pathway_by_name(pathway_name)
        if not pathways:
            print(f"未找到匹配的信号通路: {pathway_name}")
            return

        print(f"找到以下匹配的信号通路:")
        for pathway_id, name in pathways:
            print(f"  - {pathway_id}: {name}")

        if output_dir:
            output_dir = str(CommonUtils.ensure_output_dir(output_dir))
            json_output_file = Path(output_dir) / f"pathway_search_result.json"
            pathways_data = [{"pathway_id": p[0], "pathway_name": p[1]} for p in pathways]
            self._save_json_with_csv(pathways_data, json_output_file)

    def process_command(self, args: argparse.Namespace) -> None:
        """处理KEGG相关命令

        Args:
            args: 命令行参数
        """
        if args.command == "gene-pathways":
            self.process_gene_pathways(
                gene_name=args.gene_name,
                output_dir=args.output_dir
            )
        elif args.command == "download-pathway":
            self.process_download_pathway(
                pathway_id=args.pathway_id,
                output_dir=args.output_dir
            )
        elif args.command == "search-pathway":
            self.process_search_pathway(
                pathway_name=args.pathway_name,
                output_dir=args.output_dir
            )