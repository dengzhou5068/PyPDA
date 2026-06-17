# PyPDA: 蛋白质分析综合工具

PyPDA是一个集成了蛋白质序列分析、PDB文件处理和UniProt数据检索功能的综合工具，旨在为生物信息学研究提供便捷的蛋白质数据分析解决方案。

## 🚀 功能特点

- **🎯 序列分析**：从UniProt批量获取蛋白质序列、提取特定区域、执行突变分析和序列比对
- **🏗️ PDB处理**：下载PDB文件、提取配体信息并分类管理
- **📋 结构信息提取**：自动从PDB文件中提取PDB ID、结构标题、实验方法和分辨率范围，生成结构化报告
- **🧪 分子结构相似性**：基于SMILES计算小分子与PDB结构中配体的相似性，辅助分子对接
- **📊 UniProt数据检索**：获取蛋白质详细注释信息并生成结构化报告
- **🎯 OpenTargets数据查询**：查询疾病药物关联、靶点药物关联、基因疾病关联和疾病靶点关联信息
- **🧬 KEGG信号通路分析**：查询基因参与的信号通路、下载通路图、搜索通路，支持JSON和CSV双格式输出
- **🔍 增强的基因搜索功能**：实现从精确到宽松的递进式查询策略，支持灵活的基因名和蛋白质名匹配，提高检索成功率
- **⚡ 并行计算支持**：采用多线程和多进程并行处理，显著提升大规模数据处理速度
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

1. 配置文件设置：
   - `exclude_residues.ini`：定义需要排除的残基
   - `config.ini`：配置文件路径

## 🎯 使用方法

### 命令行使用

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

- 实现从精确到宽松的递进式查询策略（gene\_exact → gene → 一般搜索）
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
python pypda.py pdb fetch <query> [output_dir] [--smiles <smiles>]
```

- **参数**：
  - `query`: PDB数据库搜索查询字符串（例如：BRCA1, kinase, DNA polymerase）
  - `output_dir`: 输出目录，默认值为`result/pdb_output/查询词_YYYYMMDD_HHMMSS`
  - `--smiles`: 可选，小分子SMILES字符串，用于计算与PDB结构中配体的结构相似性

**示例**：

```bash
# 基本用法：使用蛋白质名称搜索并下载PDB文件
python pypda.py pdb fetch BRCA1

# 使用功能名称搜索
python pypda.py pdb fetch kinase

# 使用更复杂的查询
python pypda.py pdb fetch "DNA polymerase"

# 高级用法：下载PDB文件并计算结构相似性
python pypda.py pdb fetch HDAC1 --smiles "CC(=O)N(C)C(=O)N1CCC(CC1)C(C)C"  # 示例SMILES字符串
```

#### 2.3 OpenTargets数据查询工具 (`opentargets`)

用于查询OpenTargets数据库，获取疾病药物关联、基因疾病关联和疾病靶点关联信息。

##### 2.3.1 查询疾病相关药物 (`disease-drugs`)

查询特定疾病的相关药物信息

```bash
python pypda.py opentargets disease-drugs [--disease-name <疾病名称>] [--disease-id <EFO ID>] [--limit <数量>] [--format <格式>]
```

- **参数**：
  - `--disease-name`: 疾病名称（例如：breast cancer）
  - `--disease-id`: 疾病的EFO ID（例如：MONDO_0007254）
  - `--limit`: 返回结果数量限制，默认100
  - `--format`: 输出格式（json/csv/all），默认json

**示例**：

```bash
python pypda.py opentargets disease-drugs --disease-name "lung cancer" --limit 50 --format csv
```

##### 2.3.2 查询靶点相关药物 (`target-drugs`)

查询特定靶点的相关药物信息

```bash
python pypda.py opentargets target-drugs <gene_name> [--format <格式>]
```

- **参数**：
  - `gene_name`: 基因名称（例如：EGFR）
  - `--format`: 输出格式（json/csv/all），默认json

**示例**：

```bash
python pypda.py opentargets target-drugs EGFR --format csv
```

##### 2.3.3 查询基因关联疾病 (`target-associations`)

查询特定基因关联的疾病信息

```bash
python pypda.py opentargets target-associations <gene_name> [--limit <数量>] [--format <格式>]
```

- **参数**：
  - `gene_name`: 基因名称（例如：TP53）
  - `--limit`: 返回结果数量限制，默认100
  - `--format`: 输出格式（json/csv/all），默认json

**示例**：

```bash
python pypda.py opentargets target-associations EGFR --limit 50 --format csv
```

##### 2.3.4 查询疾病关联靶点 (`disease-associations`)

查询特定疾病关联的靶点信息

```bash
python pypda.py opentargets disease-associations [--disease-name <疾病名称>] [--disease-id <EFO ID>] [--limit <数量>] [--format <格式>]
```

- **参数**：
  - `--disease-name`: 疾病名称（例如：diabetes）
  - `--disease-id`: 疾病的EFO ID（例如：EFO_0010164）
  - `--limit`: 返回结果数量限制，默认100
  - `--format`: 输出格式（json/csv/all），默认json

**示例**：

```bash
python pypda.py opentargets disease-associations --disease-name "breast cancer" --limit 50 --format csv
```

#### 2.4 KEGG信号通路分析工具 (`kegg`)

用于查询KEGG数据库，获取基因参与的信号通路信息并下载通路图。

##### 2.4.1 查询基因参与的信号通路 (`gene-pathways`)

查询特定基因参与的所有信号通路，并下载对应的通路图

```bash
python pypda.py kegg gene-pathways <gene_name> [output_dir]
```

- **参数**：
  - `gene_name`: 基因名称（例如：TP53）
  - `output_dir`: 输出目录，默认值为`result/kegg_output/基因名_YYYYMMDD_HHMMSS`

**输出文件**：
- `基因名_pathways.json`：信号通路列表（JSON格式）
- `基因名_pathways.csv`：信号通路列表（CSV格式，Excel兼容）
- `pathway_images/`：信号通路图目录（PNG格式）

**示例**：

```bash
python pypda.py kegg gene-pathways TP53
```

##### 2.4.2 下载信号通路图 (`download-pathway`)

根据信号通路ID下载通路图

```bash
python pypda.py kegg download-pathway <pathway_id> [output_dir]
```

- **参数**：
  - `pathway_id`: 信号通路ID（例如：hsa04110）
  - `output_dir`: 输出目录，默认值为`result/kegg_output/通路ID_YYYYMMDD_HHMMSS`

**示例**：

```bash
python pypda.py kegg download-pathway hsa04110
```

##### 2.4.3 搜索信号通路 (`search-pathway`)

根据名称搜索信号通路

```bash
python pypda.py kegg search-pathway <pathway_name> [output_dir]
```

- **参数**：
  - `pathway_name`: 信号通路名称关键词（例如：cancer）
  - `output_dir`: 输出目录（可选）

**示例**：

```bash
python pypda.py kegg search-pathway cancer
```

#### 2.5 UniProt数据处理工具 (`uniprot`)

用于从UniProt数据库获取蛋白质注释信息并生成分析报告。

##### 2.5.1 获取UniProt数据并生成报告 (`fetch`)

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

##### 2.5.2 分析现有UniProt数据文件 (`analyze`)

分析本地已有的UniProt蛋白质信息JSON文件，提取关键数据并生成Markdown格式的分析报告。

```bash
python pypda.py uniprot analyze <file>
```

- **参数**：
  - `file`: 要分析的UniProt蛋白质信息JSON文件路径

**示例**：

```bash
python pypda.py uniprot analyze result/uniprot_reports/TP53_20240101_120000/TP53_20240101_120000.json
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
├── opentargets/           # OpenTargets查询结果
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
├── opentargets/            # OpenTargets数据查询模块
│   ├── __init__.py
│   ├── opentargets_api.py  # OpenTargets API交互
│   └── opentargets_processor.py # OpenTargets命令处理
├── kegg/                   # KEGG信号通路分析模块
│   ├── __init__.py
│   ├── kegg_api.py         # KEGG API交互
│   └── kegg_processor.py   # KEGG命令处理
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
├── pypda.py                # 主应用入口和命令行工具入口
├── requirements.txt        # 依赖列表
└── setup.py                # 安装配置
```

## ⚙️ 配置文件说明

1. **exclude_residues.ini**：定义需要排除的残基类型
   - 在配体提取时用于过滤不需要考虑的残基

## 🔧 开发信息

- **主要更新**：
  - **版本 0.6.0**：
    - 新增KEGG信号通路分析功能，支持查询基因参与的信号通路、下载通路图、搜索通路
    - KEGG模块支持JSON和CSV双格式输出，方便数据处理和分析
    - 新增OpenTargets数据查询功能，支持疾病药物查询、靶点药物查询、基因关联疾病查询和疾病关联靶点查询
    - 修复OpenTargets API查询字段问题，将`knownDrugs`更新为`drugAndClinicalCandidates`
    - 增强种属信息提取功能，支持从多个来源提取种属信息
    - 重构PDB搜索功能，改用pypdb库直接搜索PDB数据库
    - 优化PDB文件分类逻辑：先根据是否有配体分类到 no_ligand 或 with_ligand 文件夹，再根据种属进行子分类
    - 添加种属信息提取功能，将种属信息补充到 structure_info.md 文件中
    - 优化文件夹创建逻辑，避免创建空文件夹

## 📋 依赖库

核心依赖：

- biopython：生物信息学工具包
- requests：HTTP请求库
- pandas：数据处理
- configparser：配置文件解析
- pypdb：PDB数据库查询工具
- rdkit：用于分子结构处理和相似性计算（需额外安装）

完整依赖列表参见`requirements.txt`

> **注意**：RDKit库需要额外安装，可通过`pip install rdkit>=2023.03.01`命令安装。该库用于支持`--smiles`参数的结构相似性计算功能。

