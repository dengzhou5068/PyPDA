        
# Pypda - 蛋白质分析综合工具

Pypda是一个用于蛋白质序列分析和PDB文件处理的综合工具，提供序列获取、结构域分析、突变模拟、序列比对以及PDB文件下载和配体提取等功能。

## 功能特点
- 蛋白质序列获取与结构域分析
- 序列提取与突变模拟
- 多序列比对与同源性分析
- PDB文件批量下载与处理
- 配体信息提取与分析
- 并行处理优化，提高效率
- 完善的错误处理和日志记录

## 安装说明

### 依赖项
- Python 3.7+
- Biopython
- requests
- configparser

### 安装命令
```bash
pip install biopython requests configparser
```

## 配置文件
程序需要以下配置文件（位于项目根目录）：

1. **seq_config.ini** - 序列分析相关配置
2. **pdb_config.ini** - PDB处理相关配置
3. **exclude_residues.ini** - 排除残基配置

## 使用方法

### 命令行格式
```bash
python pypda.py <tool> [subcommand] [options]
```

### 可用工具

#### 1. 序列分析工具 (`seq`)
提供蛋白质序列相关的各种操作

##### 1.1 获取序列和结构域信息 (`fetch`)
从UniProt批量获取蛋白质序列和结构域信息

```bash
python pypda.py seq fetch <input_file> <output_dir>
```
- `input_file`: 包含HGNC基因名称的文本文件
- `output_dir`: 输出目录，用于保存FASTA文件和结构域信息

##### 1.2 提取子序列 (`extract`)
从FASTA文件中提取指定位置的子序列

```bash
python pypda.py seq extract <fasta_file> <start> <end>
```
- `fasta_file`: 输入FASTA文件路径
- `start`: 起始位置（1-based）
- `end`: 结束位置（1-based）

##### 1.3 序列突变 (`mutate`)
对蛋白质序列执行定点突变

```bash
python pypda.py seq mutate <fasta_file> --pos <positions> --aa <amino_acids>
```
- `fasta_file`: 输入FASTA文件路径
- `--pos`: 突变位置列表（空格分隔）
- `--aa`: 对应位置的新氨基酸（空格分隔）

##### 1.4 序列比对 (`align`)
比较多个蛋白质序列的同源性

```bash
python pypda.py seq align <fasta_files...>
```
- `fasta_files`: 多个FASTA文件路径（空格分隔）

#### 2. PDB文件处理工具 (`pdb`)
处理PDB文件下载、配体提取等任务

```bash
python pypda.py pdb
```
- 该命令会读取`pdb_config.ini`中的配置，批量处理PDB相关任务

## 输出文件说明
- 序列分析结果保存在指定的输出目录中
- PDB处理结果包括：
  - 下载的CIF文件
  - 配体信息文件（pdb_ligand.md）
  - 化学信息文件（chemical_components_info.md）
  - 配体坐标文件
  - 分类后的文件（with_ligands/和no_ligand/目录）

## 错误处理
程序错误信息会记录在`error.txt`文件中，位于各输出目录下

## 示例

### 1. 获取蛋白质序列和结构域
```bash
python pypda.py seq fetch genes.txt protein_sequences
```

### 2. 提取序列的特定区域
```bash
python pypda.py seq extract BRD4.fasta 58 457
```

### 3. 执行序列突变
```bash
python pypda.py seq mutate BRD4.fasta --pos 10 25 30 --aa A V L
```

### 4. 比较多个序列
```bash
python pypda.py seq align BRD4.fasta CDK9.fasta PCSK9.fasta
```

### 5. 处理PDB文件
```bash
python pypda.py pdb
```
        