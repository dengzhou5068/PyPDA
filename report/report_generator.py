#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
报告生成模块
"""
from pathlib import Path
from typing import Union, Dict, Any, List

from logger.logger import Logger
from uniprot.protein_analyzer import ProteinInfo


class ReportGenerator:
    """报告生成器类，负责将蛋白质信息生成为Markdown报告"""

    @staticmethod
    def generate_md_report(info: ProteinInfo, md_filename: Union[str, Path]) -> None:
        """将分析结果生成Markdown文件

        Args:
            info: 包含蛋白质信息的结构化字典
            md_filename: 要生成的Markdown文件名
        """
        try:
            with open(md_filename, 'w', encoding='utf-8') as f:
                f.write(f"# UniProt 蛋白质信息分析报告 - {info['basic_info'].get('primaryAccession', 'N/A')}\n\n")

                # 1. 基本识别信息
                f.write("## 1. 基本识别信息\n")
                bi = info['basic_info']
                f.write(f"- **条目类型**: {bi.get('entryType', 'N/A')}\n")
                f.write(f"- **主要登录号**: {bi.get('primaryAccession', 'N/A')}\n")
                if bi.get('secondaryAccessions'):
                    f.write(f"- **次要登录号**: {', '.join(bi['secondaryAccessions'])}\n")
                f.write(f"- **UniProtKB ID**: {bi.get('uniProtkbId', 'N/A')}\n")
                f.write(f"- **注释评分**: {bi.get('annotationScore', 'N/A')}\n")
                audit = bi.get('entryAudit', {})
                f.write(f"- **首次公开日期**: {audit.get('firstPublicDate', 'N/A')}\n")
                f.write(f"- **最后注释更新日期**: {audit.get('lastAnnotationUpdateDate', 'N/A')}\n")
                f.write(f"- **最后序列更新日期**: {audit.get('lastSequenceUpdateDate', 'N/A')}\n")
                f.write(f"- **条目版本**: {audit.get('entryVersion', 'N/A')}\n")
                f.write(f"- **序列版本**: {audit.get('sequenceVersion', 'N/A')}\n\n")

                # 2. 生物学背景信息
                f.write("## 2. 生物学背景信息\n")
                bio = info['biology_info']
                f.write(f"- **科学名称**: {bio.get('scientificName', 'N/A')}\n")
                f.write(f"- **常用名称**: {bio.get('commonName', 'N/A')}\n")
                f.write(f"- **分类ID**: {bio.get('taxonId', 'N/A')}\n")
                if bio.get('lineage'):
                    f.write(f"- **生物学谱系**: {' -> '.join(bio['lineage'])}\n")
                f.write(f"- **蛋白质存在证据**: {bio.get('proteinExistence', 'N/A')}\n\n")

                # 3. 蛋白质描述和基因信息
                f.write("## 3. 蛋白质描述和基因信息\n")
                pd = info['protein_desc']
                rn = pd['recommendedName']
                f.write(f"- **推荐全名**: {rn.get('fullName', 'N/A')}\n")
                if rn.get('shortNames'):
                    f.write(f"- **推荐简称**: {', '.join(rn['shortNames'])}\n")
                if rn.get('ecNumbers'):
                    f.write(f"- **EC 编号**: {', '.join(rn['ecNumbers'])}\n")

                alternative_names = pd.get('alternativeNames', [])
                if alternative_names:
                    f.write("- **备选名称**:\n")
                    for alt_name in alternative_names:
                        alt_full_name = alt_name.get('fullName', {}).get('value', 'N/A')
                        alt_ec_numbers = [ec.get('value') for ec in alt_name.get('ecNumbers', []) if ec.get('value')]
                        f.write(f"  - {alt_full_name}" + (f" (EC: {', '.join(alt_ec_numbers)})" if alt_ec_numbers else "") + "\n")

                genes = pd.get('genes', [])
                if genes:
                    for gene in genes:
                        gene_name = gene.get('geneName', {}).get('value', 'N/A')
                        f.write(f"- **基因名称**: {gene_name}\n")
                        gene_synonyms = [syn.get('value') for syn in gene.get('synonyms', []) if syn.get('value')]
                        if gene_synonyms:
                            f.write(f"  - **基因同义词**: {', '.join(gene_synonyms)}\n")
                f.write("\n")

                # 4. 功能和活性注释
                f.write("## 4. 功能和活性注释\n")
                comments_data = info['comments']
                comment_order = ['FUNCTION', 'CATALYTIC ACTIVITY', 'COFACTOR', 'ACTIVITY REGULATION',
                                 'TISSUE SPECIFICITY', 'SUBCELLULAR LOCATION', 'PTM', 'SIMILARITY',
                                 'DISEASE',]
                
                found_any_comment = False
                for ct in comment_order:
                    if comments_data.get(ct):
                        found_any_comment = True
                        f.write(f"### {ct.replace('_', ' ').title()}\n")
                        for comment in comments_data[ct]:
                            texts = []
                            if 'texts' in comment and isinstance(comment['texts'], list):
                                for text_item in comment['texts']:
                                    if isinstance(text_item, dict) and 'value' in text_item:
                                        texts.append(text_item['value'])
                                    elif isinstance(text_item, str):
                                        texts.append(text_item)
                            elif 'text' in comment:
                                if isinstance(comment['text'], dict) and 'value' in comment['text']:
                                    texts.append(comment['text']['value'])
                                elif isinstance(comment['text'], str):
                                    texts.append(comment['text'])

                            if texts:
                                f.write(f"- {'; '.join(texts)}\n")
                            
                            if ct == 'CATALYTIC ACTIVITY':
                                reaction = comment.get('reaction')
                                if reaction:
                                    f.write(f"  - **Reaction**: {reaction.get('name', 'N/A')}\n")
                            elif ct == 'DISEASE':
                                disease_name = comment.get('disease', {}).get('diseaseName', 'N/A')
                                acronym = comment.get('disease', {}).get('acronym', 'N/A')
                                f.write(f"  - **Disease**: {disease_name} ({acronym})\n")
                            elif ct == 'SUBCELLULAR LOCATION':
                                for location in comment.get('subcellularLocations', []):
                                    f.write(f"  - **Location**: {location.get('location', {}).get('value', 'N/A')}\n")
                                    for topology in location.get('topologies', []):
                                        f.write(f"    - **Topology**: {topology.get('value', 'N/A')}\n")
                                    for orientation in location.get('orientations', []):
                                        f.write(f"    - **Orientation**: {orientation.get('value', 'N/A')}\n")
                            elif ct == 'COFACTOR':
                                cofactors = comment.get('cofactors', [])
                                if cofactors:
                                    f.write("  - **Cofactors**:\n")
                                    for cofactor in cofactors:
                                        cofactor_name = cofactor.get('name', 'N/A')
                                        chebi_id = cofactor.get('cofactorCrossReference', {}).get('id', 'N/A')
                                        if chebi_id != 'N/A':
                                            f.write(f"    - {cofactor_name} (ChEBI ID: {chebi_id})\n")
                                        else:
                                            f.write(f"    - {cofactor_name}\n")
                                else:
                                    f.write("  - 未找到详细辅因子信息。\n")
           
                if not found_any_comment:
                    f.write("- 未找到功能和活性注释信息。\n")
                f.write("\n")

                # 5. 蛋白质特征
                f.write("## 5. 蛋白质特征\n")
                features = info['features']
                f.write("### 特征类型统计\n")
                if features['counts']:
                    for ft, count in features['counts'].items():
                        f.write(f"- {ft}: {count}\n")
                else:
                    f.write("- 未找到特征类型统计信息。\n")

                f.write("\n### 详细特征信息\n")
                found_any_detailed_feature = False
                detailed_feature_order = ["Chain", "Region", "Active site", "Binding site", "Modified residue", "Mutagenesis", "Domain"]
                for dt in detailed_feature_order:
                    if features['detailed'].get(dt):
                        found_any_detailed_feature = True
                        f.write(f"#### {dt}\n")
                        for feature in features['detailed'][dt]:
                            desc = feature.get('description', 'N/A')
                            if isinstance(desc, dict):
                                desc = desc.get('value', 'N/A')

                            loc = feature.get('location', {})
                            begin = loc.get('start', {}).get('value', 'N/A')
                            end = loc.get('end', {}).get('value', 'N/A')
                            
                            extra_info = []
                            if dt == 'Mutagenesis':
                                alt_seq_data = feature.get('alternativeSequence', {})
                                original_seq = alt_seq_data.get('originalSequence')
                                alternative_seqs = alt_seq_data.get('alternativeSequences', [])
                                
                                if original_seq:
                                    mutation_str = f"原始: {original_seq}"
                                    if alternative_seqs:
                                        mutation_str += f" -> 突变: {', '.join(alternative_seqs)}"
                                    extra_info.append(mutation_str)
                            
                            feature_line = f"- {desc} (位置: {begin}-{end})"
                            if extra_info:
                                feature_line += f", {', '.join(extra_info)}" 

                            f.write(feature_line + "\n")
                if not found_any_detailed_feature:
                    f.write("- 未找到详细特征信息。\n")
                f.write("\n")

                # 6. 蛋白质相互作用
                f.write("## 6. 蛋白质相互作用\n")
                interactions = info['interactions']
                if interactions:
                    for i, interaction in enumerate(interactions, 1):
                        interactant_one = interaction.get('interactantOne', {})
                        interactant_two = interaction.get('interactantTwo', {})
                        main_protein_id = interactant_one.get('uniProtKBAccession', 'N/A')
                        partner_id = interactant_two.get('uniProtKBAccession', 'N/A')
                        partner_name = interactant_two.get('geneName', interactant_two.get('name', 'N/A'))
                        organism_differ = "是" if interaction.get('organismDiffer', False) else "否"
                        f.write(f"- 相互作用 {i}: {main_protein_id} <-> **{partner_name}** (UniProt ID: {partner_id})\n")
                else:
                    f.write("- 未找到相互作用信息\n")
                f.write("\n")

                # 7. 关键词
                f.write("## 7. 关键词\n")
                keywords = info['keywords']
                if keywords:
                    f.write(f"- {', '.join(sorted(keywords))}\n")
                else:
                    f.write("- 未找到关键词信息\n")
                f.write("\n")

                # 8. 参考文献
                f.write("## 8. 参考文献\n")
                references = info['references']
                if references:
                    for i, ref in enumerate(references[:10], 1):
                        citation = ref.get('citation', {})
                        authors_list = citation.get('authors', [])
                        authors = ", ".join(authors_list) if authors_list else 'N/A'
                        title = citation.get('title', 'N/A')
                        journal = citation.get('journal', 'N/A')
                        year = citation.get('publicationDate', 'N/A')[:4]
                        pubmed_id = ""
                        for xref in citation.get('uniProtKBCrossReferences', []):
                            if xref.get('database') == 'PubMed':
                                pubmed_id = f" [PMID:{xref.get('id', '')}]"
                                break
                        f.write(f"- [{i}] {authors}, \"{title}\", *{journal}*, {year}{pubmed_id}\n")
                    if len(references) > 10:
                        f.write(f"- 显示前10篇，共{len(references)}篇参考文献\n")
                else:
                    f.write("- 未找到参考文献信息\n")
                f.write("\n")

                # 9. 交叉引用
                f.write("## 9. 交叉引用\n")
                cross_refs = info['cross_references']
                found_any_xref = False
                for db in sorted(cross_refs.keys()):
                    if cross_refs[db]:
                        found_any_xref = True
                        f.write(f"### {db}\n")
                        ids = [xref.get('id') for xref in cross_refs[db] if xref.get('id')]
                        display_ids = ids[:10]
                        f.write(f"- {', '.join(display_ids)}")
                        if len(ids) > 10:
                            f.write(f" ... (共{len(ids)}个条目)")
                        f.write("\n")
                if not found_any_xref:
                    f.write("- 未找到交叉引用信息。\n")
                f.write("\n")

                # 10. 蛋白质序列信息
                f.write("## 10. 蛋白质序列信息\n")
                sequence = info['sequence']
                f.write(f"- **序列长度**: {sequence.get('length', 'N/A')} 个氨基酸\n")
                f.write(f"- **序列版本**: {sequence.get('sequenceVersion', 'N/A')}\n")
                f.write(f"- **序列MD5**: {sequence.get('md5', 'N/A')}\n")
                seq = sequence.get('value', '')
                if seq:
                    formatted_seq = '\n'.join([seq[i:i+80] for i in range(0, len(seq), 80)])
                    f.write("- **氨基酸序列**:\n\n```\n{}\n```\n".format(formatted_seq))
                else:
                    f.write("- 未找到氨基酸序列信息。\n")

            print(f"Markdown报告已生成至 {md_filename}")
        except IOError as e:
            Logger.log_error(f"生成Markdown报告失败: {md_filename} - {e}", "error.txt")
            raise
