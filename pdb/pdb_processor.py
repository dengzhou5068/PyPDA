#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PDB文件处理模块
"""
import configparser
import json
import shutil
import sys
from functools import lru_cache
from pathlib import Path
from typing import List, Dict, Set, Any, Union, Tuple

from Bio.PDB.PDBList import PDBList
from tqdm import tqdm
import pypdb

from config.config_manager import ConfigManager
from logger.logger import Logger
from utils.common_utils import CommonUtils


try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, DataStructs
except ImportError:
    pass


# 模块级别的辅助函数，用于并行处理

def _process_single_cif_file(args: Tuple[Path, List[str]]) -> Tuple[str, List[str]]:
    """处理单个CIF文件，提取配体信息

    Args:
        args: 包含(file_path, exclude_residues)的元组

    Returns:
        (pdb_id, ligands_in_pdb)的元组
    """
    file_path, exclude_residues = args
    import warnings
    from Bio import BiopythonWarning
    from Bio.PDB.MMCIFParser import MMCIFParser
    
    # 抑制Biopython的PDBConstructionWarning警告
    warnings.filterwarnings('ignore', category=BiopythonWarning)
    
    parser = MMCIFParser()
    pdb_id = file_path.stem
    ligands_in_pdb: List[str] = []
    try:
        structure = parser.get_structure(pdb_id, str(file_path))
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.id[0].startswith('H_') and residue.resname.strip().upper() not in exclude_residues:
                        ligands_in_pdb.append(residue.resname.strip().upper())
        return pdb_id, ligands_in_pdb
    except Exception as e:
        error_msg = f"解析 {file_path.name} 时出错: {e}"
        # 注意：在进程池中无法直接访问logger，这里使用简单的print
        print(error_msg)
        return pdb_id, []


def _process_single_structure_file(file_path: Path) -> Dict[str, Any]:
    """处理单个结构文件，提取结构信息

    Args:
        file_path: 文件路径

    Returns:
        结构信息字典
    """
    pdb_id = file_path.stem
    structure_info: Dict[str, Any] = {
        'pdb_id': pdb_id,
        'title': '',
        'method': '',
        'resolution_low': '',
        'resolution_high': '',
        'species': ''
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
        
        # 提取种属信息
        import re
        
        # 1. 首先尝试从_entity_src_nat.pdbx_organism_scientific非loop结构中提取
        species_pattern = r'_entity_src_nat\.pdbx_organism_scientific\s+([\'"])([^\'"]+)\1'
        species_match = re.search(species_pattern, content, re.IGNORECASE)
        if species_match:
            structure_info['species'] = species_match.group(2)

        # 2. 如果没有找到，尝试从_entity_src_gen.pdbx_gene_src_scientific_name非loop结构中提取
        if not structure_info.get('species'):
            species_pattern = r'_entity_src_gen\.pdbx_gene_src_scientific_name\s+([\'"])([^\'"]+)\1'
            species_match = re.search(species_pattern, content, re.IGNORECASE)
            if species_match:
                structure_info['species'] = species_match.group(2)

        # 3. 如果没有找到，尝试从_pdbx_entity_src_syn.organism_scientific中提取（包括loop结构）
        if not structure_info.get('species'):
            # 尝试匹配非loop结构
            species_pattern = r'_pdbx_entity_src_syn\.organism_scientific[\s\n]+([\'"])([^\'"]+)\1'
            species_match = re.search(species_pattern, content, re.IGNORECASE)
            if species_match:
                structure_info['species'] = species_match.group(2)
            else:
                # 尝试匹配loop结构
                loop_pattern = (
                    r'loop_[\s\S]*?_pdbx_entity_src_syn\.organism_scientific'
                    r'[\s\S]*?\n([\s\S]*?)#'
                )
                entity_src_syn_match = re.search(loop_pattern, content, re.IGNORECASE)
                if entity_src_syn_match:
                    table_content = entity_src_syn_match.group(1)
                    # 在loop内容中查找物种名称
                    species_match = re.search(r'\'([A-Za-z][a-z]+\s+[A-Za-z][a-z]+)\'', table_content)
                    if species_match:
                        structure_info['species'] = species_match.group(1)

        # 4. 如果仍然没有找到，尝试从_entity_src_gen loop结构中提取
        if not structure_info.get('species'):
            loop_pattern = r'loop_[\s\S]*?_entity_src_gen\.[\s\S]*?\n([\s\S]*?)#'
            entity_src_gen_match = re.search(loop_pattern, content, re.IGNORECASE)
            if entity_src_gen_match:
                table_content = entity_src_gen_match.group(1)
                # 检查是否包含Homo sapiens、9606或human
                if re.search(r'(Homo\s+sapiens|9606|human)', table_content, re.IGNORECASE):
                    structure_info['species'] = 'Homo sapiens'
                else:
                    # 尝试从_entity_src_gen.pdbx_gene_src_scientific_name字段提取
                    scientific_name_match = re.search(r'\'([A-Za-z][a-z]+\s+[A-Za-z][a-z]+)\'', table_content)
                    if scientific_name_match:
                        structure_info['species'] = scientific_name_match.group(1)
                    # 尝试从_struct_ref.db_code字段提取（如CDK6_HUMAN）
                    if not structure_info.get('species'):
                        struct_ref_pattern = r'_struct_ref\.db_code\s+([A-Za-z0-9_]+)'
                        struct_ref_match = re.search(struct_ref_pattern, content, re.IGNORECASE)
                        if struct_ref_match:
                            db_code = struct_ref_match.group(1)
                            # 检查是否包含_HUMAN后缀
                            if '_HUMAN' in db_code:
                                structure_info['species'] = 'Homo sapiens'
        
        return structure_info
    except Exception as e:
        error_msg = f"解析 {file_path.name} 时出错: {e}"
        # 注意：在进程池中无法直接访问logger，这里使用简单的print
        print(error_msg)
        return structure_info


def _process_single_pocket_file(args: Tuple[Path, List[str]]) -> List[Tuple[str, str, str]]:
    """处理单个文件，分析口袋信息

    Args:
        args: 包含(file_path, exclude_residues)的元组

    Returns:
        口袋分析结果列表
    """
    file_path, exclude_residues = args
    import warnings
    from Bio import BiopythonWarning
    from Bio.PDB import PDBParser
    from Bio.PDB.MMCIFParser import MMCIFParser
    from Bio.PDB.Selection import unfold_entities
    from Bio.PDB.NeighborSearch import NeighborSearch
    
    # 抑制Biopython的PDBConstructionWarning警告
    warnings.filterwarnings('ignore', category=BiopythonWarning)
    
    file_name = file_path.name
    pdb_id = file_path.stem
    results: List[Tuple[str, str, str]] = []
    
    try:
        # 根据文件扩展名选择解析器
        if file_path.suffix.lower() == '.cif':
            parser = MMCIFParser()
        else:  # .pdb
            parser = PDBParser()
        
        structure = parser.get_structure(pdb_id, str(file_path))
        
        # 获取所有原子用于邻居搜索
        all_atoms = unfold_entities(structure, 'A')
        neighbor_search = NeighborSearch(all_atoms)
        
        # 遍历结构中的配体
        for model in structure:
            for chain in model:
                for residue in chain:
                    # 检查是否为配体（HETATM且不在排除列表中）
                    if residue.id[0].startswith('H_') and residue.resname.strip().upper() not in exclude_residues:
                        ligand_name = residue.resname.strip().upper()
                        
                        # 获取配体的所有原子
                        ligand_atoms = list(residue.get_atoms())
                        
                        # 收集配体周围4.5埃内的残基
                        pocket_residues: Set[str] = set()
                        
                        for atom in ligand_atoms:
                            # 搜索4.5埃内的所有原子
                            nearby_atoms = neighbor_search.search(atom.coord, 4.5, level='A')
                            
                            # 过滤出氨基酸残基（非HETATM且不在排除列表中）
                            for nearby_atom in nearby_atoms:
                                nearby_residue = nearby_atom.get_parent()
                                # 检查是否为氨基酸残基（非HETATM且不在排除列表中）
                                if not nearby_residue.id[0].startswith('H_'):
                                    res_name = nearby_residue.resname.strip().upper()
                                    # 排除非氨基酸残基（如溶剂分子）
                                    if res_name not in exclude_residues:
                                        # 格式化残基为XXXyyy格式
                                        res_num = nearby_residue.id[1]
                                        formatted_residue = f"{res_name}{res_num}"
                                        pocket_residues.add(formatted_residue)
                        
                        # 将残基列表排序
                        if pocket_residues:
                            sorted_residues = sorted(list(pocket_residues))
                            residues_str = ', '.join(sorted_residues)
                            results.append((file_name, ligand_name, residues_str))
        return results
    except Exception as e:
        error_msg = f"分析 {file_name} 时出错: {e}"
        # 注意：在进程池中无法直接访问logger，这里使用简单的print
        print(error_msg)
        return [(file_name, "错误", str(e))]


class PDBProcessor:
    """PDB处理类，负责PDB文件相关操作"""
    def __init__(self, config: ConfigManager, logger: Logger, uniprot_api: Any):
        """初始化PDBProcessor实例
        
        Args:
            config: 配置管理器实例
            logger: 日志记录器实例
            uniprot_api: UniProtAPI实例，用于访问UniProt数据
        """
        self.config = config
        self.logger = logger
        self.uniprot_api = uniprot_api
        self.exclude_residues = self.parse_exclude_residues()
        self.rdkit_available = 'rdkit' in sys.modules

    def parse_exclude_residues(self) -> List[str]:
        """解析排除残基的配置文件 `exclude_residues.ini`
        
        期望配置文件中有一个或多个 [residues] 或类似 section，
        每个 key 对应一个或多个逗号分隔的残基名称。

        Returns:
            排除的残基名称列表，全部转换为大写
        """
        config = configparser.ConfigParser()
        # 使用绝对路径引用 exclude_residues.ini 文件
        config_path = Path(__file__).parent.parent / 'exclude_residues.ini'
        if not config_path.exists():
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

    def get_pdb_ids_from_query(self, query: str, max_results: int = 5) -> List[str]:
        """通过查询字符串搜索PDB ID列表
        使用 pypdb 库搜索 PDB 数据库。

        Args:
            query: 搜索查询字符串
            max_results: 最大返回结果数量，默认 5 个

        Returns:
            PDB ID列表
        """
        try:
            # 使用 pypdb 搜索 PDB ID
            pdb_ids = pypdb.Query(query).search()
            if pdb_ids:
                unique_ids = list(set(pdb_ids))
                return unique_ids[:max_results]
            else:
                self.logger.log_error(f"未找到与查询 '{query}' 匹配的PDB ID。", self.config.error_log_path)
                return []
        except Exception as e:
            self.logger.log_error(f"通过查询 '{query}' 搜索PDB ID时出错: {e}", self.config.error_log_path)
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
        
        ligand_pdb_dict: Dict[str, List[str]] = {}

        folder_path = Path(folder_path)
        if not folder_path.is_dir():
            self.logger.log_error(f"指定的PDB文件夹路径不存在或不是目录: {folder_path}", self.config.error_log_path)
            return {}

        cif_files = list(folder_path.glob('*.cif'))
        if not cif_files:
            self.logger.log_error(f"在 {folder_path} 中未找到任何.cif文件进行配体提取。", self.config.error_log_path)
            return {}

        # 准备参数列表
        args_list = [(file_path, self.exclude_residues) for file_path in cif_files]

        # 并行处理所有文件
        results = []
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor() as executor:
            for result in executor.map(_process_single_cif_file, args_list):
                results.append(result)

        # 合并结果
        for pdb_id, ligands in results:
            for ligand in set(ligands):
                ligand_pdb_dict.setdefault(ligand, []).append(pdb_id)

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
        """根据配体信息移动CIF文件到相应的子文件夹

        Args:
            folder_path: 包含CIF文件的文件夹路径
            ligand_pdb_dict: 配体到PDB ID的映射字典
        """
        folder_path = Path(folder_path)
        no_ligand_dir = folder_path / 'no_ligand'
        with_ligands_dir = folder_path / 'with_ligand'
        
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
                target_dir_name = 'with_ligand'
            else:
                dst_path = no_ligand_dir / src_path.name
                target_dir_name = 'no_ligand'
            
            try:
                if src_path != dst_path:
                    shutil.move(src_path, dst_path)
            except Exception as e:
                self.logger.log_error(f"移动文件 {src_path.name} 到 {target_dir_name} 时出错: {e}", self.config.error_log_path)

    def move_files_based_on_species(self, folder_path: Union[str, Path]) -> None:
        """根据种属信息对CIF文件进行分类

        Args:
            folder_path: 包含CIF文件的文件夹路径
        """
        folder_path = Path(folder_path)
        
        # 检查 no_ligand 和 with_ligand 文件夹
        no_ligand_dir = folder_path / 'no_ligand'
        with_ligands_dir = folder_path / 'with_ligand'
        
        # 搜索所有CIF文件
        cif_files = []
        if no_ligand_dir.exists():
            cif_files.extend(list(no_ligand_dir.glob('*.cif')))
        if with_ligands_dir.exists():
            cif_files.extend(list(with_ligands_dir.glob('*.cif')))
        
        if not cif_files:
            self.logger.log_error(f"在 {folder_path} 中未找到任何.cif文件进行种属分类。", self.config.error_log_path)
            return
        
        # 提取种属信息
        structure_info_list = self.extract_structure_info(folder_path)
        species_pdb_dict = {}
        for info in structure_info_list:
            species = info.get('species', 'Unknown')
            pdb_id = info.get('pdb_id')
            if pdb_id:
                species_pdb_dict.setdefault(species, []).append(pdb_id)
        
        # 为每个种属创建文件夹并移动文件
        for species, pdb_ids in species_pdb_dict.items():
            # 清理种属名称，去除特殊字符
            safe_species_name = (
                species.replace(' ', '_').replace('/', '_').replace('\\', '_')
                .replace(':', '_').replace('*', '_').replace('?', '_')
                .replace('"', '_').replace('<', '_').replace('>', '_')
                .replace('|', '_')
            )
            
            # 跳过空种属名称
            if not species:
                continue
                
            for pdb_id in pdb_ids:
                # 查找文件位置
                for cif_file in cif_files:
                    if cif_file.stem == pdb_id:
                        src_path = cif_file
                        # 确定目标目录
                        if src_path.parent == no_ligand_dir:
                            if no_ligand_dir.exists():
                                species_dir = no_ligand_dir / f'species_{safe_species_name}'
                                # 只在需要时创建文件夹
                                species_dir.mkdir(exist_ok=True)
                                dst_path = species_dir / src_path.name
                        elif src_path.parent == with_ligands_dir:
                            if with_ligands_dir.exists():
                                species_dir = with_ligands_dir / f'species_{safe_species_name}'
                                # 只在需要时创建文件夹
                                species_dir.mkdir(exist_ok=True)
                                dst_path = species_dir / src_path.name
                        else:
                            continue
                        
                        try:
                            if src_path != dst_path:
                                shutil.move(src_path, dst_path)
                        except Exception as e:
                            log_msg = (
                                f"移动文件 {src_path.name} 到 {species_dir.name} 时出错: {e}"
                            )
                            self.logger.log_error(log_msg, self.config.error_log_path)
                        break
        
        print(f"已根据种属对 {len(cif_files)} 个CIF文件进行分类")

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

        CommonUtils.parallel_executor(
            download_single_ligand_json,
            list(unique_ligands),
            description="Downloading ligand JSONs"
        )

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
                        log_msg = f"JSON file for {ligand} not found in {json_dir}"
                        self.logger.log_error(log_msg, self.config.error_log_path)
                    except Exception as e:
                        log_msg = f"读取或解析配体 {ligand} 的JSON文件失败: {e}"
                        self.logger.log_error(log_msg, self.config.error_log_path)
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
            log_msg = f"在 {folder_path} 中未找到任何.cif文件进行结构信息提取。"
            self.logger.log_error(log_msg, self.config.error_log_path)
            return structure_info_list
        
        # 并行处理所有文件
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor() as executor:
            results = list(executor.map(_process_single_structure_file, cif_files))

        # 合并结果
        structure_info_list.extend(results)
        
        return structure_info_list
    
    def write_structure_info_to_md(
        self,
        structure_info_list: List[Dict[str, Any]],
        output_dir: Union[str, Path]
    ) -> None:
        """将结构信息写入structure_info.md文件

        Args:
            structure_info_list: 结构信息列表
            output_dir: 输出目录
        """
        output_file_path = Path(output_dir) / 'structure_info.md'
        try:
            with open(output_file_path, 'w', encoding='utf-8') as outfile:
                outfile.write('# 结构信息汇总\n\n')
                outfile.write('| PDB ID | 结构标题 | 实验方法 | 分辨率范围 (Å) | 种属 |\n')
                outfile.write('| --- | --- | --- | --- | --- |\n')
                
                for structure_info in sorted(
                    structure_info_list,
                    key=lambda x: x['pdb_id']
                ):
                    pdb_id = structure_info['pdb_id']
                    title = structure_info['title']
                    method = structure_info['method']
                    resolution_low = structure_info['resolution_low']
                    resolution_high = structure_info['resolution_high']
                    species = structure_info.get('species', '')
                    
                    if resolution_low and resolution_high:
                        resolution = f"{resolution_low} - {resolution_high}"
                    elif resolution_low:
                        resolution = resolution_low
                    elif resolution_high:
                        resolution = resolution_high
                    else:
                        resolution = ""
                    
                    outfile.write(f'| {pdb_id} | {title} | {method} | {resolution} | {species} |\n')
            
            print(f'结构信息已写入 {output_file_path}')
        except IOError as e:
            log_msg = (
                f"写入结构信息到MD文件失败: {output_file_path} - {e}"
            )
            self.logger.log_error(log_msg, self.config.error_log_path)
            raise

    def calculate_similarity_and_rank_pdbs(
        self,
        user_smiles: str,
        output_dir: str,
        top_n: int = 5
    ) -> List[Tuple[str, float]]:
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

    def process(
        self,
        query: str = None,
        uniprot_id: str = None,
        output_dir: str = None,
        user_smiles: str = None
    ) -> None:
        """处理PDB文件下载和分析

        Args:
            query: 搜索查询字符串（优先使用）
            uniprot_id: UniProt蛋白质编号（备选）
            output_dir: 输出目录路径
            user_smiles: 可选，用户输入的小分子SMILES字符串，用于计算结构相似性
        """
        output_dir_path = Path(output_dir)
        output_dir_path.mkdir(parents=True, exist_ok=True)

        if not query and not uniprot_id:
            self.logger.log_error("缺少查询字符串或 'uniprot' ID。", self.config.error_log_path)
            return
        
        pdb_ids = []
        if query:
            print(f"\n--- 开始处理 PDB 任务 (查询: {query}) ---")
            pdb_ids = self.get_pdb_ids_from_query(query)
            if pdb_ids:
                print(f"为查询 '{query}' 找到 PDB IDs: {', '.join(pdb_ids)}")
                self.download_pdb_files(pdb_ids, output_dir)
            else:
                self.logger.log_error(f"未找到与查询 '{query}' 匹配的PDB ID。", self.config.error_log_path)
                return
        else:
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
            log_msg = (
                f"开始计算分子结构相似性，用户输入的SMILES: {user_smiles}"
            )
            self.logger.log_info(log_msg, self.config.info_log_path)
            self.calculate_similarity_and_rank_pdbs(user_smiles, str(output_dir))
        
        # 提取结构信息并写入structure_info.md文件
        structure_info_list = self.extract_structure_info(output_dir)
        if structure_info_list:
            self.write_structure_info_to_md(structure_info_list, output_dir)

        # 根据种属进行分类
        self.move_files_based_on_species(output_dir)

        print('\nPDB任务处理完成。')

    def analyze_pockets(self, folder_path: str, output_dir: str) -> None:
        """对指定文件夹下的PDB或CIF文件进行口袋分析，列出配体周围4.5埃内的氨基酸残基

        Args:
            folder_path: 包含PDB或CIF文件的文件夹路径
            output_dir: 输出目录，用于保存分析结果
        """
        import warnings
        from Bio import BiopythonWarning
        
        # 抑制Biopython的PDBConstructionWarning警告
        warnings.filterwarnings('ignore', category=BiopythonWarning)
        
        folder_path = Path(folder_path)
        output_dir_path = Path(output_dir)
        output_dir_path.mkdir(parents=True, exist_ok=True)
        
        # 查找所有PDB和CIF文件
        pdb_files = list(folder_path.glob('*.pdb')) + list(folder_path.glob('*.cif'))
        
        if not pdb_files:
            print(f"在 {folder_path} 中未找到任何PDB或CIF文件。")
            return
        
        print(f"找到 {len(pdb_files)} 个PDB/CIF文件，开始口袋分析...")
        
        # 准备参数列表
        args_list = [(file_path, self.exclude_residues) for file_path in pdb_files]

        # 并行处理所有文件
        all_results = []
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor() as executor:
            for results in executor.map(_process_single_pocket_file, args_list):
                all_results.extend(results)

        # 准备Markdown输出文件
        output_md_file = output_dir_path / 'pocket_analysis.md'
        
        with open(output_md_file, 'w', encoding='utf-8') as md_file:
            md_file.write('# 口袋分析结果\n\n')
            md_file.write('| 文件名 | 配体 | 口袋残基 (XXXyyy格式) |\n')
            md_file.write('| --- | --- | --- |\n')
            
            for result in all_results:
                file_name, ligand_name, residues_str = result
                md_file.write(f"| {file_name} | {ligand_name} | {residues_str} |\n")
        
        print(f"口袋分析完成，结果已写入 {output_md_file}")
