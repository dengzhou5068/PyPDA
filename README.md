# PyPDA: 蛋白质分析综合工具

PyPDA是一个集成了蛋白质序列分析、PDB文件处理和UniProt数据检索功能的综合工具，旨在为生物信息学研究提供便捷的蛋白质数据分析解决方案。

## 🚀 功能特点

- **🎯 序列分析**：从UniProt批量获取蛋白质序列、提取特定区域、执行突变分析和序列比对
- **🏗️ PDB处理**：下载PDB文件、提取配体信息并分类管理
- **📋 结构信息提取**：自动从PDB文件中提取PDB ID、结构标题、实验方法和分辨率范围，生成结构化报告
- **🧪 分子结构相似性**：基于SMILES计算小分子与PDB结构中配体的相似性，辅助分子对接
- **📊 UniProt数据检索**：获取蛋白质详细注释信息并生成结构化报告
- **🔍 增强的基因搜索功能**：实现从精确到宽松的递进式查询策略，支持灵活的基因名和蛋白质名匹配，提高检索成功率
- **🌐 Gradio界面**：提供直观的Web界面，无需命令行即可使用所有功能
- **📁 统一管理**：所有结果文件和下载zip文件统一存储在`result/`目录下

## 📦 安装方法

### 前提条件
- Python 3.8+
- 所需依赖库：参见requirements.txt

### 安装步骤

1. 克隆或下载本项目到本地

2. 安装依赖包：
```bash
pip install -r requirements.txt
```

3. 配置文件设置：
   - `exclude_residues.ini`：定义需要排除的残基
   - `config.ini`：配置文件路径

## 🎯 使用方法

### 1. Gradio Web界面（推荐）

启动直观的Web界面：
```bash
python launch_gradio.py
```

访问 http://localhost:7860 即可使用所有功能。

### 2. 命令行使用

#### 基本命令格式
```bash
python pypda.py <工具> <命令> [参数]
```

#### 2.1 序列分析工具 (`seq`)

用于蛋白质序列的获取、提取、突变和比对分析。

##### 2.1.1 获取蛋白质序列和结构域信息 (`fetch`)
从UniProt批量获取蛋白质序列和结构域信息，支持增强的搜索策略，提高检索成功率

```bash
python pypda.py seq fetch <genes> [output_dir]
```
- **参数**：
  - `genes`: 要下载的基因名称，用空格分隔（例如：BRCA1 TP53 EGFR）
  - `output_dir`: 输出目录，默认值为`result/protein_sequences/基因名_YYYYMMDD_HHMMSS`

**增强的搜索功能**：
- 实现从精确到宽松的递进式查询策略（gene_exact → gene → 一般搜索）
- 支持灵活的基因名和蛋白质名匹配（精确匹配和包含匹配）
- 多结果筛选，提高匹配准确度
- 完善的错误处理和日志记录

**示例**：
```bash
python pypda.py seq fetch BRCA1 TP53 EGFR
```

即使是非标准或难以匹配的基因名（如MOGT2、MOGT3）也能成功检索

##### 2.1.2 提取序列区域 (`extract`)
从FASTA文件中提取指定位置的氨基酸序列

```bash
python pypda.py seq extract <fasta_file> <start> <end> [--output_dir <目录>]
```
- **参数**：
  - `fasta_file`: 输入FASTA文件路径
  - `start`: 起始位置（1-based）
  - `end`: 结束位置（1-based）
  - `--output_dir`: 输出目录，默认值为输入文件所在目录

**示例**：
```bash
python pypda.py seq extract protein.fasta 10 50 --output_dir ./extract_result
```

##### 2.1.3 执行序列突变 (`mutate`)
对蛋白质序列执行定点突变

```bash
python pypda.py seq mutate <fasta_file> --pos <位置列表> --aa <氨基酸列表> [--output_dir <目录>]
```
- **参数**：
  - `fasta_file`: 输入FASTA文件路径
  - `--pos`: 突变位置列表（空格分隔）
  - `--aa`: 对应位置的新氨基酸（空格分隔）
  - `--output_dir`: 输出目录，默认值为输入文件所在目录

**示例**：
```bash
python pypda.py seq mutate protein.fasta --pos 15 23 --aa A K --output_dir ./mutate_result
```

##### 2.1.4 序列比对分析 (`align`)
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

#### 2.2 PDB文件处理工具 (`pdb`)

用于PDB文件的下载、配体提取和分类管理，以及基于小分子SMILES的结构相似性计算。

```bash
python pypda.py pdb fetch <accession> [output_dir] [--smiles <smiles>]
```
- **参数**：
  - `accession`: 蛋白质名称或基因名称（例如：HDAC1）
  - `output_dir`: 输出目录，默认值为`result/pdb_output/基因名_YYYYMMDD_HHMMSS`
  - `--smiles`: 可选，小分子SMILES字符串，用于计算与PDB结构中配体的结构相似性

**示例**：
```bash
# 基本用法：下载PDB文件并提取配体信息
python pypda.py pdb fetch HDAC1

# 高级用法：下载PDB文件并计算结构相似性
python pypda.py pdb fetch HDAC1 --smiles "CC(=O)N(C)C(=O)N1CCC(CC1)C(C)C"  # 示例SMILES字符串
```

#### 2.3 UniProt数据处理工具 (`uniprot`)

用于从UniProt数据库获取蛋白质注释信息并生成分析报告。

##### 2.3.1 获取UniProt数据并生成报告 (`fetch`)
从UniProt API获取蛋白质数据并生成JSON和Markdown报告

```bash
python pypda.py uniprot fetch <accession> [output_dir]
```
- **参数**：
  - `accession`: 蛋白质名称或基因名称（例如：TP53）
  - `output_dir`: 输出目录，默认值为`result/uniprot_reports/基因名_YYYYMMDD_HHMMSS`

**示例**：
```bash
python pypda.py uniprot fetch TP53
```

## 📁 目录结构说明

项目采用统一的目录结构管理所有输出文件：

```
result/
├── zip/                    # 所有下载zip文件统一存放
├── extract_output/         # 序列提取结果
├── mutate_output/          # 序列突变结果
├── align_output/           # 序列比对结果
├── pdb_output/             # PDB处理结果
├── protein_sequences/      # 蛋白质序列结果
├── uniprot_reports/        # UniProt数据结果
└── temp/                   # 临时文件目录
```

## 📊 输出文件说明

### 序列分析输出
- **FASTA文件**：包含蛋白质序列
- **domain_info.md**：结构域信息报告
- **JSON文件**：原始UniProt数据（命名格式：`基因名_YYYYMMDD_HHMMSS.json`）
- **Markdown报告**：蛋白质详细信息（命名格式：`基因名_YYYYMMDD_HHMMSS.md`）

### PDB处理输出
- **PDB文件**：蛋白质结构文件
- **配体文件**：提取的配体信息
- **分类报告**：配体分类结果
- **结构信息报告**：`structure_info.md`，包含PDB ID、结构标题、实验方法和分辨率范围
- **结构相似性报告**：当使用`--smile`参数时生成，包含与输入小分子最相似的PDB结构排名

## 📁 项目结构

PyPDA采用模块化设计，代码结构清晰，便于维护和扩展：

```
pypda/
├── config/                 # 配置管理模块
│   ├── __init__.py
│   └── config_manager.py   # 配置文件管理
├── logger/                 # 日志管理模块
│   ├── __init__.py
│   └── logger.py           # 日志记录功能
├── pdb/                    # PDB处理模块
│   ├── __init__.py
│   └── pdb_processor.py    # PDB文件处理功能
├── report/                 # 报告生成模块
│   ├── __init__.py
│   └── report_generator.py # 报告生成功能
├── sequence/               # 序列处理模块
│   ├── __init__.py
│   └── sequence_processor.py # 序列分析功能
├── uniprot/                # UniProt数据模块
│   ├── __init__.py
│   ├── protein_analyzer.py # 蛋白质数据分析
│   └── uniprot_api.py      # UniProt API交互
├── utils/                  # 工具函数模块
│   ├── __init__.py
│   └── common_utils.py     # 通用工具函数
├── gradio_app.py           # Gradio Web界面
├── launch_gradio.py        # Gradio启动脚本
├── pypda.py                # 主应用入口和命令行工具入口
├── requirements.txt        # 依赖列表
└── setup.py                # 安装配置
```

## ⚙️ 配置文件说明

1. **exclude_residues.ini**：定义需要排除的残基类型
   - 在配体提取时用于过滤不需要考虑的残基

## 🔧 开发信息

- **主要文件**：
  - `pypda.py`：主应用入口和命令行工具入口（原 main.py 已合并至此）
  - `gradio_app.py`：Gradio Web界面
  - `launch_gradio.py`：Gradio启动脚本
  - `config/config_manager.py`：配置管理
  - `uniprot/uniprot_api.py`：UniProt API交互
  - `sequence/sequence_processor.py`：序列分析功能
  - `pdb/pdb_processor.py`：PDB文件处理

- **主要更新**：
  - 优化了`search_uniprot_by_name`方法，实现递进式查询策略
  - 增强了`fetch_protein_sequences`方法的序列描述验证逻辑
  - 修复了Gradio界面与后端目录结构一致性问题
  - 修复了UniProt API字段名问题，将'name'改为'protein_name'
  - 采用模块化设计，将代码拆分为多个功能模块，提高可维护性
  - 统一了命令行参数格式，将位置参数改为可选参数（--output_dir）
  - 将`main.py`合并到`pypda.py`，简化项目结构，统一入口点
  - 添加了结构信息提取功能，自动生成`structure_info.md`文件，包含PDB ID、结构标题、实验方法和分辨率范围

## 📋 依赖库

核心依赖：
- biopython：生物信息学工具包
- requests：HTTP请求库
- gradio：Web界面框架
- pandas：数据处理
- configparser：配置文件解析
- rdkit：用于分子结构处理和相似性计算（需额外安装）

完整依赖列表参见`requirements.txt`

> **注意**：RDKit库需要额外安装，可通过`pip install rdkit>=2023.03.01`命令安装。该库用于支持`--smiles`参数的结构相似性计算功能。