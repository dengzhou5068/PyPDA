#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
蛋白质数据分析模块
"""
from typing import Dict, Any, List, TypedDict


class ProteinInfo(TypedDict):
    """蛋白质信息结构化数据类型定义

    包含从UniProt JSON数据中提取的各类蛋白质信息
    """
    basic_info: Dict[str, Any]
    biology_info: Dict[str, Any]
    protein_desc: Dict[str, Any]
    comments: Dict[str, List[Dict[str, Any]]]
    features: Dict[str, Any]
    interactions: List[Dict[str, Any]]
    keywords: List[str]
    references: List[Dict[str, Any]]
    cross_references: Dict[str, List[Dict[str, Any]]]
    sequence: Dict[str, Any]


class ProteinAnalyzer:
    """蛋白质数据分析类，负责从JSON数据中提取和组织蛋白质信息"""

    @staticmethod
    def extract_protein_info(data: Dict[str, Any]) -> ProteinInfo:
        """从JSON数据中提取蛋白质信息并组织为结构化字典

        Args:
            data: 从UniProt API获取的原始JSON数据

        Returns:
            包含各类蛋白质信息的结构化字典
        """
        info: ProteinInfo = {
            'basic_info': {},
            'biology_info': {},
            'protein_desc': {},
            'comments': {},
            'features': {},
            'interactions': [],
            'keywords': [],
            'references': [],
            'cross_references': {},
            'sequence': {}
        }

        # 1. 基本识别信息
        info['basic_info'] = {
            'entryType': data.get('entryType', 'N/A'),
            'primaryAccession': data.get('primaryAccession', 'N/A'),
            'secondaryAccessions': data.get('secondaryAccessions', []),
            'uniProtkbId': data.get('uniProtkbId', 'N/A'),
            'annotationScore': data.get('annotationScore', 'N/A'),
            'entryAudit': data.get('entryAudit', {})
        }

        # 2. 生物学背景信息
        organism = data.get('organism', {})
        info['biology_info'] = {
            'scientificName': organism.get('scientificName', 'N/A'),
            'commonName': organism.get('commonName', 'N/A'),
            'taxonId': organism.get('taxonId', 'N/A'),
            'lineage': organism.get('lineage', []),
            'proteinExistence': data.get('proteinExistence', 'N/A')
        }

        # 3. 蛋白质描述和基因信息
        protein_desc = data.get('proteinDescription', {})
        recommended_name = protein_desc.get('recommendedName', {})
        info['protein_desc'] = {
            'recommendedName': {
                'fullName': recommended_name.get('fullName', {}).get('value', 'N/A'),
                'shortNames': [sn.get('value') for sn in recommended_name.get('shortNames', []) if sn.get('value')],
                'ecNumbers': [ec.get('value') for ec in recommended_name.get('ecNumbers', []) if ec.get('value')]
            },
            'alternativeNames': protein_desc.get('alternativeNames', []),
            'genes': data.get('genes', [])
        }

        # 4. 功能和活性注释
        comments = data.get('comments', [])
        comment_types = ['FUNCTION', 'CATALYTIC ACTIVITY', 'COFACTOR', 'ACTIVITY REGULATION',
                         'TISSUE SPECIFICITY', 'SUBCELLULAR LOCATION', 'PTM', 'SIMILARITY', 'DISEASE', 'INTERACTION']
        comment_type_map = {ct.lower(): ct for ct in comment_types}
        info['comments'] = {ct: [] for ct in comment_types}

        for comment in comments:
            ct_lower = comment.get('commentType', '').lower()
            if ct_lower in comment_type_map:
                ct = comment_type_map[ct_lower]
                info['comments'][ct].append(comment)
                if ct == 'INTERACTION':
                    info['interactions'].extend(comment.get('interactions', []))

        # 5. 蛋白质特征
        features = data.get('features', [])
        info['features'] = {
            'counts': {},
            'detailed': {}
        }
        detailed_feature_types = [
            "Chain", "Region", "Active site", "Binding site",
            "Modified residue", "Mutagenesis", "Domain"
        ]
        for dt in detailed_feature_types:
            info['features']['detailed'][dt] = []

        for feature in features:
            ft = feature.get('type')
            info['features']['counts'][ft] = info['features']['counts'].get(ft, 0) + 1
            if ft == "Region" and feature.get('description') == "Domain":
                info['features']['detailed']['Domain'].append(feature)
            elif ft in info['features']['detailed']:
                info['features']['detailed'][ft].append(feature)

        # 6. 蛋白质相互作用(已在4. 功能和活性注释中提取)

        # 7. 关键词
        info['keywords'] = [kw.get('name') for kw in data.get('keywords', []) if kw.get('name')]

        # 8. 参考文献
        info['references'] = data.get('references', [])

        # 9. 交叉引用
        cross_references = data.get('uniProtKBCrossReferences', [])
        important_databases = ["PDB", "DrugBank", "GO", "Reactome", "HGNC", "GeneID", "KEGG", "AlphaFoldDB", "IntAct"]
        info['cross_references'] = {db: [] for db in important_databases}
        for xref in cross_references:
            db = xref.get('database')
            if db in info['cross_references']:
                info['cross_references'][db].append(xref)

        # 10. 蛋白质序列信息
        info['sequence'] = data.get('sequence', {})

        return info
