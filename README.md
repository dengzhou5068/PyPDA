# Python Protein Data Analysis 工具说明
## 项目概述
pypda.py是一个综合性的蛋白质分析工具，整合了序列分析和PDB结构分析功能，能够帮助研究人员快速获取蛋白质序列、分析结构域、进行序列突变和比对，以及提取PDB文件中的配体信息。

## 功能特性
### 序列分析功能
- 从UniProt数据库获取蛋白质序列
- 提取蛋白质结构域信息并保存为Markdown文件
- 截取蛋白质序列的特定片段
- 对蛋白质序列进行定点突变
- 多序列比对并计算同源性百分比
### PDB结构分析功能
- 从UniProt ID获取相关的PDB ID
- 下载PDB文件（mmCIF格式）
- 提取PDB文件中的配体信息
- 分离包含配体和不包含配体的PDB文件
- 下载配体的化学信息（JSON格式）
- 生成配体信息汇总（Markdown格式）
- 提取配体坐标并保存为单独的CIF文件
## 安装要求
- Python 3.x
- 依赖库：BioPython, requests, configparser
## 配置文件
项目需要以下配置文件：

1. seq_config.ini ：序列分析相关配置
   
   - UNIPROT_API：UniProt API地址
   - UNIPROT_XML_API：UniProt XML格式API地址
   - XML_NAMESPACES：XML命名空间配置
2. pdb_config.ini ：PDB分析相关配置
   
   - 包含任务特定的输出目录和UniProt ID
3. exclude_residues.ini ：需要排除的残基配置
## 使用方法
### 命令行格式
```
python pypda.py [seq|pdb] [command] [options]
```
### 序列分析示例
1. 获取蛋白质序列和结构域信息：
```
python pypda.py seq fetch --input_file genes.txt --output_dir seq_output
```
2. 提取序列片段（从10到100位）：
```
python pypda.py seq extract --fasta_file protein.fasta --start 10 --end 100
```
3. 序列突变（在位置5和15分别突变为A和K）：
```
python pypda.py seq mutate --fasta_file protein.fasta --pos 5 15 --aa A K
```
4. 多序列比对：
```
python pypda.py seq align --fasta_files protein1.fasta protein2.fasta
```
### PDB分析示例
```
python pypda.py pdb
```
## 输出文件说明
- 序列分析：FASTA文件、结构域信息Markdown文件
- PDB分析：
  - 配体信息：pdb_ligand.md
  - 化学组件信息：chemical_components_info.md
  - 分类文件：with_ligands/和no_ligand/目录
  - 配体坐标文件：ligands/目录下的CIF文件
## 注意事项
- 确保配置文件路径正确
- 网络连接稳定以保证成功下载数据
- 大型PDB文件可能需要较长处理时间