# PyPDA: 蛋白质分析综合工具

PyPDA是一个集成了蛋白质序列分析、PDB文件处理和UniProt数据检索功能的综合工具，旨在为生物信息学研究提供便捷的蛋白质数据分析解决方案。

## 功能特点

- **序列分析**：从UniProt批量获取蛋白质序列、提取特定区域、执行突变分析和序列比对
- **PDB处理**：下载PDB文件、提取配体信息并分类管理
- **UniProt数据检索**：获取蛋白质详细注释信息并生成结构化报告
- **多格式输出**：支持FASTA、JSON和Markdown等多种格式的输入输出

## 安装方法

### 前提条件
- Python 3.8+ 
- 所需依赖库：biopython, requests, configparser

### 安装步骤

1. 克隆或下载本项目到本地

2. 安装依赖包

3. 配置文件设置：
   - exclude_residues.ini：定义需要排除的残基

## 使用说明

### 基本命令格式
```bash
python pypda.py <工具> <命令> [参数]
```

### 工具和命令详解

#### 1. 序列分析工具 (`seq`)

用于蛋白质序列的获取、提取、突变和比对分析。

##### 1.1 获取蛋白质序列和结构域信息 (`fetch`)
从UniProt批量获取蛋白质序列和结构域信息

```bash
python pypda.py seq fetch <genes> <output_dir>
```
- **参数**：
  - `genes`: 要下载的基因名称，用空格分隔（例如：BRCA1 TP53 EGFR）
  - `output_dir`: 输出目录，用于保存FASTA序列文件和结构域信息报告，默认值为"protein_sequences"

**示例**：
```bash
python pypda.py seq fetch genes.txt ./sequences
```

##### 1.2 提取序列区域 (`extract`)
从FASTA文件中提取指定位置的氨基酸序列

```bash
python pypda.py seq extract <fasta_file> <start> <end>
```
- **参数**：
  - `fasta_file`: 输入FASTA文件路径
  - `start`: 起始位置（1-based）
  - `end`: 结束位置（1-based）

**示例**：
```bash
python pypda.py seq extract protein.fasta 10 50
```

##### 1.3 执行序列突变 (`mutate`)
对蛋白质序列执行定点突变

```bash
python pypda.py seq mutate <fasta_file> --pos <位置列表> --aa <氨基酸列表>
```
- **参数**：
  - `fasta_file`: 输入FASTA文件路径
  - `--pos`: 突变位置列表（空格分隔）
  - `--aa`: 对应位置的新氨基酸（空格分隔）

**示例**：
```bash
python pypda.py seq mutate protein.fasta --pos 15 23 --aa A K
```

##### 1.4 序列比对分析 (`align`)
比较多个蛋白质序列的同源性

```bash
python pypda.py seq align <fasta_files>
```
- **参数**：
  - `fasta_files`: 多个FASTA文件路径（空格分隔）

**示例**：
```bash
python pypda.py seq align protein1.fasta protein2.fasta protein3.fasta
```

#### 2. PDB文件处理工具 (`pdb`)

用于PDB文件的下载、配体提取和分类管理。
```bash
python pypda pdb <accession> pdb_output
```
- **参数**：
  - `accession`: 蛋白质名称或基因名称（例如：HDAC1）
  - `output_dir`: 输出目录，默认值为"pdb_output"

#### 3. UniProt数据处理工具 (`uniprot`)

用于从UniProt数据库获取蛋白质注释信息并生成分析报告。

##### 3.1 获取UniProt数据并生成报告 (`fetch`)
从UniProt API获取蛋白质数据并生成JSON和Markdown报告

```bash
python pypda.py uniprot fetch <accession>
```
- **参数**：
  - `accession`: 蛋白质名称或基因名称（例如：TP53）

**示例**：
```bash
python pypda.py uniprot fetch TP53
```

##### 3.2 分析JSON文件并生成报告 (`analyze`)
分析已有的UniProt JSON数据文件并生成Markdown报告

```bash
python pypda.py uniprot analyze -f <json_file>
```
- **参数**：
  - `-f, --file`: 要分析的JSON文件路径

**示例**：
```bash
python pypda.py uniprot analyze -f <UniProt_ID>_<date>.json
```

## 输出文件说明

- **序列分析输出**：
  - FASTA文件：包含蛋白质序列
  - domain_info.md：结构域信息报告

- **UniProt分析输出**：
  - JSON文件：原始数据（命名格式：<UniProt_ID>_<date>.json）
  - Markdown报告：蛋白质详细信息（命名格式：<UniProt_ID>_<date>.md）

## 配置文件说明

1. **exclude_residues.ini**：定义需要排除的残基类型
   - 在配体提取时用于过滤不需要考虑的残基