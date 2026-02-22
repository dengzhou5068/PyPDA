# 📁 PyPDA 项目结构与文件存储

## 🎯 项目结构说明

PyPDA 是一个蛋白质数据分析综合工具，包含命令行界面和 Gradio Web 界面。项目采用模块化设计，功能清晰，易于扩展和维护。

## 🗂️ 完整项目结构

```
pypda/
├── pypda.py                        # 主程序入口（原 main.py 已合并至此）
├── setup.py                        # 项目安装配置
├── requirements.txt                # 依赖包列表
├── README.md                       # 项目说明文档
├── FOLDER_STRUCTURE.md             # 本文件，项目结构说明
├── exclude_residues.ini            # 排除残基配置文件
├── error.txt                       # 错误日志文件
├── .gitignore                      # Git 忽略文件配置
├── __init__.py                     # 包初始化文件
├── config/                         # 配置管理模块
│   ├── config_manager.py           # 配置管理类
│   └── __init__.py
├── logger/                         # 日志管理模块
│   ├── logger.py                   # 日志记录类
│   └── __init__.py
├── sequence/                       # 序列处理模块
│   ├── sequence_processor.py       # 序列处理核心逻辑
│   └── __init__.py
├── pdb/                            # PDB 文件处理模块
│   ├── pdb_processor.py            # PDB 处理核心逻辑
│   └── __init__.py
├── uniprot/                        # UniProt 数据处理模块
│   ├── uniprot_api.py              # UniProt API 调用
│   ├── protein_analyzer.py         # 蛋白质数据分析
│   └── __init__.py
├── report/                         # 报告生成模块
│   ├── report_generator.py         # 报告生成逻辑
│   └── __init__.py
├── utils/                          # 通用工具模块
│   ├── common_utils.py             # 通用工具函数
│   └── __init__.py
└── result/                         # 所有结果文件存储目录
    ├── protein_sequences/          # 蛋白质序列获取结果
    │   └── protein_YYYYMMDD_HHMMSS/ # 按时间戳命名的批次文件夹
    ├── extract_output/             # 序列提取结果
    │   └── extract_YYYYMMDD_HHMMSS/ # 按时间戳命名的批次文件夹
    ├── mutate_output/              # 序列突变结果
    │   └── mutate_YYYYMMDD_HHMMSS/ # 按时间戳命名的批次文件夹
    ├── align_output/               # 序列比对结果
    │   └── align_YYYYMMDD_HHMMSS/  # 按时间戳命名的批次文件夹
    ├── uniprot_reports/            # UniProt 数据获取结果
    │   └── uniprot_YYYYMMDD_HHMMSS/ # 按时间戳命名的批次文件夹
    ├── pdb_output/                 # PDB 文件处理结果
    │   └── pdb_YYYYMMDD_HHMMSS/    # 按时间戳命名的批次文件夹
    │       ├── chemical_components_info.md # 化学组分信息报告
    │       ├── pdb_ligand.md       # PDB 配体信息报告
    │       ├── json/               # 配体 JSON 数据
    │       ├── no_ligand/          # 无配体的 PDB 文件
    │       └── with_ligands/       # 有配体的 PDB 文件
    └── zip/                        # 压缩文件存储
        └── *.zip                   # 结果压缩包
```

## ✨ 存储系统特点

### 1. **按功能分类存储**
- 每个功能模块的输出文件存储在独立的文件夹中
- 清晰区分不同功能的结果，便于查找和管理

### 2. **时间戳批次管理**
- 每次操作自动创建带时间戳的子文件夹
- 格式：`功能名_YYYYMMDD_HHMMSS`
- 支持重复操作，不会覆盖之前的结果

### 3. **PDB 结果详细分类**
- PDB 文件按是否含配体分类存储
- 配体信息单独存储为 JSON 格式
- 生成详细的配体和化学组分报告

### 4. **统一结果管理**
- 所有结果文件集中在 `result/` 目录下
- 支持结果压缩，便于分享和保存

## 🚀 使用方法

### 命令行方式
```bash
python pypda.py --help             # 查看帮助信息
python pypda.py seq fetch BRCA1    # 获取 BRCA1 蛋白质序列
python pypda.py pdb TP53           # 处理 TP53 的 PDB 文件
python pypda.py uniprot fetch EGFR # 获取 EGFR 的 UniProt 数据
```

## 📍 文件访问

- **项目根目录**：`pypda\`
- **结果文件**：`pypda\result\`

## 📌 最近更新

- **主入口合并**：将 `main.py` 合并到 `pypda.py`，简化项目结构
- **存储优化**：新增 `zip` 文件夹用于存储压缩结果
- **PDB 处理增强**：按配体情况分类存储 PDB 文件

现在您可以轻松管理和查找不同批次的结果文件了！🎉