# Pypda: 蛋白质分析综合工具

Pypda是一个功能强大的蛋白质分析工具，集成了序列处理和PDB文件分析功能，帮助研究人员高效处理蛋白质数据。

## 功能特点

- **序列处理**：批量获取蛋白质序列、提取结构域、执行序列突变、多序列比对
- **PDB分析**：获取PDB ID、下载PDB文件、提取配体信息、处理配体数据
- **统一配置管理**：通过配置文件管理API和路径信息
- **并行处理**：支持多线程下载和处理，提高效率

## 安装说明

### 前提条件

- Python 3.7+
- 依赖库：biopython, requests

### 安装步骤

1. 克隆或下载项目到本地
2. 安装依赖：
   ```bash
   pip install biopython requests
   ```

## 使用方法

### 命令行接口

Pypda提供两个主要工具：`seq`（序列处理）和`pdb`（PDB分析）。

#### 序列处理工具 (`seq`)

```bash
python pypda.py seq <command> [options]
```

可用命令：

1. **fetch**: 批量获取蛋白质序列和结构域信息
   ```bash
   python pypda.py seq fetch <input_file> <output_dir>
   ```
   - `input_file`: 包含HGNC基因名称的文本文件（每行一个基因）
   - `output_dir`: 输出目录，默认值为"protein_sequences"

2. **extract**: 提取序列的指定区域
   ```bash
   python pypda.py seq extract <fasta_file> <start> <end>
   ```
   - `fasta_file`: 输入FASTA文件路径
   - `start`: 起始位置（1-based）
   - `end`: 结束位置（1-based）

3. **mutate**: 执行蛋白质序列突变
   ```bash
   python pypda.py seq mutate <fasta_file> --pos <positions> --aa <amino_acids>
   ```
   - `fasta_file`: 输入FASTA文件路径
   - `--pos`: 突变位置列表（空格分隔）
   - `--aa`: 新氨基酸列表（空格分隔）

4. **align**: 比较多个序列
   ```bash
   python pypda.py seq align <fasta_files>
   ```
   - `fasta_files`: 多个FASTA文件路径（空格分隔）

#### PDB分析工具 (`pdb`)

目前PDB工具正在开发中，敬请期待。

## 配置文件

项目使用以下配置文件：

- `seq_config.ini`: 序列处理相关配置
- `pdb_config.ini`: PDB分析相关配置
- `exclude_residues.ini`: 排除残基配置

## 示例

### 获取蛋白质序列和结构域信息

```bash
# 创建包含基因名称的文件
echo -e "BRD4\nPCSK9" > genes.txt

# 批量获取序列和结构域信息
python pypda.py seq fetch genes.txt brd4_pcsk9_sequences
```

### 提取序列区域

```bash
python pypda.py seq extract brd4_pcsk9_sequences/BRD4.fasta 58 457
```

### 执行序列突变

```bash
python pypda.py seq mutate BRD4.fasta --pos 10 20 --aa A F
```

### 多序列比对

```bash
python pypda.py seq align BRD4.fasta PCSK9.fasta
```

## 依赖项

- biopython: 用于生物信息学分析
- requests: 用于HTTP请求
- configparser: 用于配置文件解析