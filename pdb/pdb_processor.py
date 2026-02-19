#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PDB文件处理模块
"""
import configparser
import json
import os
import shutil
import sys
from functools import lru_cache
from pathlib import Path
from typing import List, Dict, Set, Any, Union, Tuple

from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBList import PDBList
from tqdm import tqdm

from config.config_manager import ConfigManager
from logger.logger import Logger
from utils.common_utils import CommonUtils


try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, DataStructs
except ImportError:
    pass


class PDBProcessor:
    """PDB处理类，负责PDB文件相关操作"""
    def __init__(self, config: ConfigManager, logger: Logger, uniprot_api: Any):
        self.config = config
        self.logger = logger
        self.uniprot_api = uniprot_api
        self.exclude_residues = self.parse_exclude_residues()
        self.rdkit_available = 'rdkit' in sys.modules

    def parse_exclude_residues(self) -> List[str]:
        """解析排除残基的配置文件 `exclude_residues.ini`。
        期望配置文件中有一个或多个 [residues] 或类似 section，
        每个 key 对应一个或多个逗号分隔的残基名称。
        """
        config = configparser.ConfigParser()
        config_path = 'exclude_residues.ini'
        if not Path(config_path).exists():
            self.logger.log_error(f"配置文件 {config_path} 未找到。将不排除任何残基。", self.config.error_log_path)
            return []
        
        try:
            config.read(config_path)
        except Exception as e:
            self.logger.log_error(f"读取 {config_path} 配置文件失败: {e}", self.config.error_log_path)
            return []

        exclude_residues = []
        for section in config.sections():
            for key in config[section]:
                exclude_residues.extend([res.strip().upper() for res in config[section][key].split(',') if res.strip()])
        return list(set(exclude_residues))

    @lru_cache(maxsize=None)
    def get_pdb_ids_from_uniprot(self, uniprot_id: str) -> List[str]:
        """通过UniProt ID获取对应的PDB ID列表（带缓存优化）
        使用 UniProt REST API 获取交叉引用信息。

        Args:
            uniprot_id: UniProt ID

        Returns:
            PDB ID列表
        """
        try:
            data = self.uniprot_api.get_uniprot_data(uniprot_id)
            if not data:
                self.logger.log_error(f"无法从UniProt API获取UniProt ID {uniprot_id} 的数据。", self.config.error_log_path)
                return []
            
            pdb_ids = []
            cross_references = data.get('uniProtKBCrossReferences', [])
            for xref in cross_references:
                if xref.get('database') == 'PDB':
                    pdb_ids.append(xref.get('id'))
            return list(set(pdb_ids))
        except Exception as e:
            self.logger.log_error(f"获取UniProt ID {uniprot_id} 的PDB ID时出错: {e}", self.config.error_log_path)
            return []

    def download_pdb_files(self, pdb_ids: List[str], output_dir: str) -> None:
        """下载PDB文件

        Args:
            pdb_ids: PDB ID列表
            output_dir: 输出目录
        """
        output_dir_path = Path(output_dir)
        output_dir_path.mkdir(parents=True, exist_ok=True)
        pdbl = PDBList()
        failed_files: List[str] = []

        def download_single_file(pdb_id: str) -> None:
            try:
                local_path = pdbl.retrieve_pdb_file(pdb_id, pdir=str(output_dir_path), file_format='mmCif')
                if local_path and Path(local_path).exists():
                    pass
                else:
                    self.logger.log_error(f"下载CIF文件 {pdb_id} 失败，未返回有效路径或文件不存在。", self.config.error_log_path)
                    failed_files.append(pdb_id)
            except Exception as e:
                error_msg = f"下载CIF文件 {pdb_id} 时出错: {e}"
                self.logger.log_error(error_msg, self.config.error_log_path)
                failed_files.append(pdb_id)

        CommonUtils.parallel_executor(download_single_file, pdb_ids, description="Downloading PDB files")

        if failed_files:
            failed_log_path = output_dir_path / 'download_failed.txt'
            with open(failed_log_path, 'a', encoding='utf-8') as f:
                for failed_file_id in failed_files:
                    f.write(f"{failed_file_id}.cif\n")
            print(f"部分CIF文件下载失败，详情请查看 {failed_log_path}")

    def extract_ligands_from_pdb(self, folder_path: Union[str, Path]) -> Dict[str, List[str]]:
        """从.cif文件中提取配体信息

        Args:
            folder_path: 文件夹路径

        Returns:
            配体到PDB ID的映射字典
        """
        import warnings
        from Bio import BiopythonWarning
        # 抑制Biopython的PDBConstructionWarning警告
        warnings.filterwarnings('ignore', category=BiopythonWarning)
        
        parser = MMCIFParser()
        ligand_pdb_dict: Dict[str, List[str]] = {}

        folder_path = Path(folder_path)
        if not folder_path.is_dir():
            self.logger.log_error(f"指定的PDB文件夹路径不存在或不是目录: {folder_path}", self.config.error_log_path)
            return {}

        cif_files = list(folder_path.glob('*.cif'))
        if not cif_files:
            self.logger.log_error(f"在 {folder_path} 中未找到任何.cif文件进行配体提取。", self.config.error_log_path)
            return {}

        for file_path in tqdm(cif_files, desc="Extracting ligands"):
            pdb_id = file_path.stem
            try:
                structure = parser.get_structure(pdb_id, str(file_path))
                ligands_in_pdb: List[str] = []
                for model in structure:
                    for chain in model:
                        for residue in chain:
                            if residue.id[0].startswith('H_') and residue.resname.strip().upper() not in self.exclude_residues:
                                ligands_in_pdb.append(residue.resname.strip().upper())
                
                for ligand in set(ligands_in_pdb):
                    ligand_pdb_dict.setdefault(ligand, []).append(pdb_id)

            except Exception as e:
                error_msg = f"解析 {file_path.name} 时出错: {e}"
                self.logger.log_error(error_msg, self.config.error_log_path)
                continue
        return ligand_pdb_dict

    @staticmethod
    def write_ligand_info_to_md(ligand_pdb_dict: Dict[str, List[str]], output_dir: Union[str, Path]) -> None:
        """将配体信息写入pdb_ligand.md文件

        Args:
            ligand_pdb_dict: 配体到PDB ID的映射字典
            output_dir: 输出目录
        """
        output_file_path = Path(output_dir) / 'pdb_ligand.md'
        try:
            with open(output_file_path, 'w', encoding='utf-8') as outfile:
                outfile.write('| Ligands | PDB ID |\n')
                outfile.write('| --- | --- |\n')
                for ligand, pdb_ids in sorted(ligand_pdb_dict.items()):
                    pdb_id_str = ', '.join(sorted(list(set(pdb_ids))))
                    outfile.write(f'| {ligand} | {pdb_id_str} |\n')
            print(f'结果已写入 {output_file_path}')
        except IOError as e:
            Logger.log_error(f"写入配体信息到MD文件失败: {output_file_path} - {e}", "error.txt")
            raise

    def move_files_based_on_ligands(self, folder_path: Union[str, Path], ligand_pdb_dict: Dict[str, List[str]]) -> None:
        """根据配体信息移动CIF文件到'no_ligand'或'with_ligands'子文件夹

        Args:
            folder_path: 包含CIF文件的文件夹路径
            ligand_pdb_dict: 配体到PDB ID的映射字典
        """
        folder_path = Path(folder_path)
        no_ligand_dir = folder_path / 'no_ligand'
        with_ligands_dir = folder_path / 'with_ligands'

        no_ligand_dir.mkdir(exist_ok=True)
        with_ligands_dir.mkdir(exist_ok=True)

        pdb_ids_with_ligands = set()
        for pdb_ids in ligand_pdb_dict.values():
            pdb_ids_with_ligands.update(pdb_ids)

        cif_files = list(folder_path.glob('*.cif'))
        for src_path in tqdm(cif_files, desc="Moving PDB files"):
            pdb_id = src_path.stem
            
            if pdb_id in pdb_ids_with_ligands:
                dst_path = with_ligands_dir / src_path.name
                target_dir_name = 'with_ligands'
            else:
                dst_path = no_ligand_dir / src_path.name
                target_dir_name = 'no_ligand'
            
            try:
                if src_path != dst_path:
                    shutil.move(src_path, dst_path)
            except Exception as e:
                self.logger.log_error(f"移动文件 {src_path.name} 到 {target_dir_name} 时出错: {e}", self.config.error_log_path)

    def download_ligand_json(self, unique_ligands: Set[str], output_dir: Union[str, Path]) -> None:
        """通过API并行查询配体的JSON文件并下载保存到json子文件夹

        Args:
            unique_ligands: 唯一配体集合
            output_dir: 输出目录
        """
        json_dir = Path(output_dir) / 'json'
        json_dir.mkdir(exist_ok=True)

        def download_single_ligand_json(ligand: str) -> None:
            session = self.uniprot_api.session
            try:
                url = f"https://data.rcsb.org/rest/v1/core/chemcomp/{ligand}"
                response = session.get(url, timeout=10)
                response.raise_for_status()
                
                data = response.json()
                json_path = json_dir / f'{ligand}.json'
                with open(json_path, 'w', encoding='utf-8') as json_file:
                    json.dump(data, json_file, ensure_ascii=False, indent=4)
            except Exception as e:
                if hasattr(e, 'response') and e.response.status_code == 404:
                    self.logger.log_error(f"配体 {ligand} 未在RCSB PDB找到 (404 Not Found)。", self.config.error_log_path)
                else:
                    self.logger.log_error(f"下载配体 {ligand} 时出错: {e}", self.config.error_log_path)

        CommonUtils.parallel_executor(download_single_ligand_json, list(unique_ligands), description="Downloading ligand JSONs")

    def write_chemical_info_to_md(self, unique_ligands: Set[str], output_dir: Union[str, Path]) -> None:
        """从json子文件夹中读取JSON文件并写入chemical_components_info.md

        Args:
            unique_ligands: 唯一配体集合
            output_dir: 输出目录
        """
        json_dir = Path(output_dir) / 'json'
        output_md_file = Path(output_dir) / 'chemical_components_info.md'
        
        try:
            with open(output_md_file, 'w', encoding='utf-8') as md_file:
                md_file.write("| Chemical Component ID | name | formula | formula_weight | Canonical Smiles |\n")
                md_file.write("| --- | --- | --- | --- | --- |\n")
                for ligand in tqdm(sorted(list(unique_ligands)), desc="Writing chemical info"):
                    json_path = json_dir / f'{ligand}.json'
                    try:
                        data = CommonUtils.load_from_json(str(json_path))
                        chem_comp = data.get('chem_comp', {})
                        name = chem_comp.get('name', '')
                        formula = chem_comp.get('formula', '')
                        formula_weight = chem_comp.get('formula_weight', '')
                        canonical_smiles = ''
                        descriptors = data.get('pdbx_chem_comp_descriptor', [])
                        for descriptor in descriptors:
                            if descriptor.get('type') == 'SMILES_CANONICAL':
                                canonical_smiles = descriptor.get('descriptor', '')
                                break
                        md_file.write(f"| {ligand} | {name} | {formula} | {formula_weight} | {canonical_smiles} |\n")
                    except FileNotFoundError:
                        self.logger.log_error(f"JSON file for {ligand} not found in {json_dir}", self.config.error_log_path)
                    except Exception as e:
                        self.logger.log_error(f"读取或解析配体 {ligand} 的JSON文件失败: {e}", self.config.error_log_path)
            print(f"MD file {output_md_file} has been created successfully.")
        except IOError as e:
            self.logger.log_error(f"创建或写入 chemical_components_info.md 文件失败: {e}", self.config.error_log_path)
            raise
    
    def extract_structure_info(self, folder_path: Union[str, Path]) -> List[Dict[str, Any]]:
        """从.cif文件中提取结构信息

        Args:
            folder_path: 包含CIF文件的文件夹路径

        Returns:
            结构信息列表，每个元素包含PDB ID、结构标题、实验方法、分辨率范围
        """
        structure_info_list = []
        folder_path = Path(folder_path)
        
        # 搜索所有子文件夹中的CIF文件
        cif_files = list(folder_path.rglob('*.cif'))
        if not cif_files:
            self.logger.log_error(f"在 {folder_path} 中未找到任何.cif文件进行结构信息提取。", self.config.error_log_path)
            return structure_info_list
        
        for file_path in tqdm(cif_files, desc="Extracting structure info"):
            pdb_id = file_path.stem
            structure_info = {
                'pdb_id': pdb_id,
                'title': '',
                'method': '',
                'resolution_low': '',
                'resolution_high': ''
            }
            
            try:
                with open(file_path, 'r', encoding='utf-8') as f:
                    content = f.read()
                
                # 提取PDB ID
                import re
                
                # 提取结构标题
                title_match = re.search(r'_struct.title\s+([^\n]+)', content)
                if title_match:
                    title = title_match.group(1).strip().strip("'\"")
                    structure_info['title'] = title
                
                # 提取实验方法
                method_match = re.search(r'_exptl.method\s+([^\n]+)', content)
                if method_match:
                    method = method_match.group(1).strip().strip("'\"")
                    structure_info['method'] = method
                
                # 提取分辨率范围 - 优先使用refine中的值，其次使用reflns中的值
                refine_res_low_match = re.search(r'_refine.ls_d_res_low\s+([^\n]+)', content)
                refine_res_high_match = re.search(r'_refine.ls_d_res_high\s+([^\n]+)', content)
                if refine_res_low_match and refine_res_high_match:
                    structure_info['resolution_low'] = refine_res_low_match.group(1).strip()
                    structure_info['resolution_high'] = refine_res_high_match.group(1).strip()
                else:
                    reflns_res_low_match = re.search(r'_reflns.d_resolution_low\s+([^\n]+)', content)
                    reflns_res_high_match = re.search(r'_reflns.d_resolution_high\s+([^\n]+)', content)
                    if reflns_res_low_match and reflns_res_high_match:
                        structure_info['resolution_low'] = reflns_res_low_match.group(1).strip()
                        structure_info['resolution_high'] = reflns_res_high_match.group(1).strip()
                
                structure_info_list.append(structure_info)
            except Exception as e:
                error_msg = f"解析 {file_path.name} 时出错: {e}"
                self.logger.log_error(error_msg, self.config.error_log_path)
                continue
        
        return structure_info_list
    
    def write_structure_info_to_md(self, structure_info_list: List[Dict[str, Any]], output_dir: Union[str, Path]) -> None:
        """将结构信息写入structure_info.md文件

        Args:
            structure_info_list: 结构信息列表
            output_dir: 输出目录
        """
        output_file_path = Path(output_dir) / 'structure_info.md'
        try:
            with open(output_file_path, 'w', encoding='utf-8') as outfile:
                outfile.write('# 结构信息汇总\n\n')
                outfile.write('| PDB ID | 结构标题 | 实验方法 | 分辨率范围 (Å) |\n')
                outfile.write('| --- | --- | --- | --- |\n')
                
                for structure_info in sorted(structure_info_list, key=lambda x: x['pdb_id']):
                    pdb_id = structure_info['pdb_id']
                    title = structure_info['title']
                    method = structure_info['method']
                    resolution_low = structure_info['resolution_low']
                    resolution_high = structure_info['resolution_high']
                    
                    if resolution_low and resolution_high:
                        resolution = f"{resolution_low} - {resolution_high}"
                    elif resolution_low:
                        resolution = resolution_low
                    elif resolution_high:
                        resolution = resolution_high
                    else:
                        resolution = ""
                    
                    outfile.write(f'| {pdb_id} | {title} | {method} | {resolution} |\n')
            
            print(f'结构信息已写入 {output_file_path}')
        except IOError as e:
            self.logger.log_error(f"写入结构信息到MD文件失败: {output_file_path} - {e}", self.config.error_log_path)
            raise

    def calculate_similarity_and_rank_pdbs(self, user_smiles: str, output_dir: str, top_n: int = 5) -> List[Tuple[str, float]]:
        """计算用户输入的小分子SMILES与PDB结构中配体的结构相似性，并排序PDB结构

        Args:
            user_smiles: 用户输入的小分子SMILES字符串
            output_dir: 输出目录路径
            top_n: 返回前N个最高相似性的PDB结构

        Returns:
            排序后的PDB ID和相似性分数的元组列表
        """
        if not self.rdkit_available:
            self.logger.log_error("RDKit库未正确安装，无法计算结构相似性。", self.config.error_log_path)
            return []

        # 读取chemical_components_info.md文件获取配体SMILES信息
        md_file_path = Path(output_dir) / 'chemical_components_info.md'
        if not md_file_path.exists():
            self.logger.log_error(f"配体信息文件 {md_file_path} 不存在，无法获取配体SMILES信息。", self.config.error_log_path)
            return []

        try:
            # 读取用户输入的分子
            user_mol = Chem.MolFromSmiles(user_smiles)
            if user_mol is None:
                self.logger.log_error(f"无效的SMILES字符串: {user_smiles}", self.config.error_log_path)
                return []
            # 使用RDKit计算分子指纹
            user_fp = AllChem.GetMorganFingerprintAsBitVect(user_mol, 2, nBits=2048)

            # 读取配体信息
            ligand_smiles_dict = {}
            with open(md_file_path, 'r', encoding='utf-8') as f:
                lines = f.readlines()
            for line in lines[2:]:  # 跳过表头和分隔线
                parts = [p.strip() for p in line.strip().split('|') if p.strip()]
                if len(parts) >= 5 and parts[4]:  # 检查是否有SMILES信息
                    ligand_id = parts[0]
                    smiles = parts[4]
                    ligand_smiles_dict[ligand_id] = smiles

            # 读取pdb_ligand.md获取配体-PDB映射
            pdb_ligand_file = Path(output_dir) / 'pdb_ligand.md'
            ligand_pdb_mapping = {}
            if pdb_ligand_file.exists():
                with open(pdb_ligand_file, 'r', encoding='utf-8') as f:
                    lines = f.readlines()
                for line in lines[2:]:  # 跳过表头和分隔线
                    parts = [p.strip() for p in line.strip().split('|') if p.strip()]
                    if len(parts) >= 2:
                        ligand_id = parts[0]
                        pdb_ids = [pdb_id.strip() for pdb_id in parts[1].split(',')]
                        ligand_pdb_mapping[ligand_id] = pdb_ids

            # 计算相似性并排序
            similarity_results = []
            processed_pdbs = set()  # 避免重复的PDB ID

            for ligand_id, pdb_ids in ligand_pdb_mapping.items():
                if ligand_id in ligand_smiles_dict:
                    ligand_smiles = ligand_smiles_dict[ligand_id]
                    try:
                        ligand_mol = Chem.MolFromSmiles(ligand_smiles)
                        if ligand_mol is not None:
                            # 计算配体分子指纹
                            ligand_fp = AllChem.GetMorganFingerprintAsBitVect(ligand_mol, 2, nBits=2048)
                            similarity = DataStructs.TanimotoSimilarity(user_fp, ligand_fp)
                            
                            for pdb_id in pdb_ids:
                                if pdb_id not in processed_pdbs:
                                    similarity_results.append((pdb_id, similarity))
                                    processed_pdbs.add(pdb_id)
                    except Exception as e:
                        self.logger.log_error(f"计算配体 {ligand_id} 相似性时出错: {e}", self.config.error_log_path)
                        continue

            # 根据相似性分数排序，取前N个
            similarity_results.sort(key=lambda x: x[1], reverse=True)
            top_results = similarity_results[:top_n]
            
            # 将结果写入文件
            result_file = Path(output_dir) / 'structure_similarity_ranking.md'
            with open(result_file, 'w', encoding='utf-8') as f:
                f.write("# 分子结构相似性排序结果\n\n")
                f.write(f"用户输入的SMILES: {user_smiles}\n\n")
                f.write("| 排名 | PDB ID | 结构相似性分数 (Tanimoto系数) |\n")
                f.write("| --- | --- | --- |\n")
                for i, (pdb_id, score) in enumerate(top_results, 1):
                    f.write(f"| {i} | {pdb_id} | {score:.4f} |\n")
            
            print(f"已将结构相似性排序结果写入 {result_file}")
            print(f"前{min(len(top_results), top_n)}个最相似的PDB结构：")
            for i, (pdb_id, score) in enumerate(top_results, 1):
                print(f"{i}. {pdb_id}: {score:.4f}")
            
            return top_results
        except Exception as e:
            self.logger.log_error(f"计算结构相似性时出错: {e}", self.config.error_log_path)
            return []

    def process(self, uniprot_id: str, output_dir: str, user_smiles: str = None) -> None:
        """处理PDB文件下载和分析

        Args:
            uniprot_id: UniProt蛋白质编号
            output_dir: 输出目录路径
            user_smiles: 可选，用户输入的小分子SMILES字符串，用于计算结构相似性
        """
        output_dir_path = Path(output_dir)
        output_dir_path.mkdir(parents=True, exist_ok=True)

        if not uniprot_id:
            self.logger.log_error("缺少 'uniprot' ID。", self.config.error_log_path)
            return
        
        print(f"\n--- 开始处理 PDB 任务 (UniProt ID: {uniprot_id}) ---")
        
        pdb_ids = self.get_pdb_ids_from_uniprot(uniprot_id)
        if pdb_ids:
            print(f"为 UniProt ID {uniprot_id} 找到 PDB IDs: {', '.join(pdb_ids)}")
            self.download_pdb_files(pdb_ids, output_dir)
        else:
            self.logger.log_error(f"未找到 UniProt ID {uniprot_id} 对应的PDB ID。", self.config.error_log_path)
            return

        downloaded_cif_files = list(Path(output_dir).glob('*.cif'))
        if not downloaded_cif_files:
            self.logger.log_error(f"在 {output_dir} 中未找到下载的CIF文件，跳过后续配体处理。", self.config.error_log_path)
            return

        ligand_pdb_dict = self.extract_ligands_from_pdb(output_dir)
        
        if ligand_pdb_dict:
            self.write_ligand_info_to_md(ligand_pdb_dict, output_dir)
        else:
            print(f"在 {output_dir} 中未找到任何配体信息。")
        
        self.move_files_based_on_ligands(output_dir, ligand_pdb_dict)
        md_file_path = Path(output_dir) / 'pdb_ligand.md'
        unique_ligands: Set[str] = set()
        if md_file_path.exists():
            try:
                with open(md_file_path, 'r', encoding='utf-8') as f:
                    lines = f.readlines()
                for line in lines[2:]:
                    parts = [p.strip() for p in line.strip().split('|') if p.strip()]
                    if len(parts) >= 1:
                        unique_ligands.add(parts[0])
            except IOError as e:
                self.logger.log_error(f"读取 {md_file_path} 失败: {e}", self.config.error_log_path)
        else:
            self.logger.log_error(f"配体信息文件 {md_file_path} 不存在，无法获取唯一配体列表。", self.config.error_log_path)

        if unique_ligands:
            print(f"找到唯一配体: {', '.join(sorted(list(unique_ligands)))}")
            self.download_ligand_json(unique_ligands, output_dir)
            self.write_chemical_info_to_md(unique_ligands, output_dir)
        else:
            print("未找到任何配体，跳过配体JSON下载和坐标提取。")

        # 如果提供了用户SMILES，计算结构相似性并排序PDB结构
        if user_smiles:
            self.logger.log_info(f"开始计算分子结构相似性，用户输入的SMILES: {user_smiles}", self.config.info_log_path)
            self.calculate_similarity_and_rank_pdbs(user_smiles, str(output_dir))
        
        # 提取结构信息并写入structure_info.md文件
        structure_info_list = self.extract_structure_info(output_dir)
        if structure_info_list:
            self.write_structure_info_to_md(structure_info_list, output_dir)

        print('\nPDB任务处理完成。')
