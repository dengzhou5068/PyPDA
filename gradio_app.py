#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PyPDA Gradio Web界面
提供用户友好的Web界面来访问PyPDA的各种功能
"""

import gradio as gr
import os
import sys
import tempfile
import shutil
import subprocess

# 将当前目录添加到Python路径
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from pypda import ConfigManager, Logger, CommonUtils

class PyPDAGradioInterface:
    """Gradio界面类"""
    
    def __init__(self):
        self.temp_dir = "result"
        self.zip_dir = os.path.join("result", "zip")
        os.makedirs(self.zip_dir, exist_ok=True)
        self.config = ConfigManager()
        self.logger = Logger()
        
    def cleanup_temp_files(self):
        """清理临时文件"""
        try:
            shutil.rmtree(self.temp_dir, ignore_errors=True)
            self.temp_dir = tempfile.mkdtemp(prefix="pypda_")
        except Exception as e:
            print(f"清理临时文件时出错: {e}")
    
    def get_protein_sequences(self, genes_text: str, progress=gr.Progress()):
        """获取蛋白质序列和结构域信息"""
        if not genes_text.strip():
            return "请输入基因名称", None, None
            
        genes = [g.strip() for g in genes_text.split() if g.strip()]
        if not genes:
            return "请输入有效的基因名称", None, None
        
        # 使用CommonUtils.get_output_dir生成输出目录和时间戳
        output_dir, timestamp = CommonUtils.get_output_dir("result", "protein_sequences", "protein")
        
        try:
            progress(0.1, desc="正在获取蛋白质序列...")
            
            # 构建命令
            cmd = [
                sys.executable, "pypda.py", "seq", "fetch"
            ] + genes + [output_dir]
            
            # 执行命令
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=os.path.dirname(os.path.abspath(__file__)))
            
            progress(0.7, desc="正在处理结果...")
            
            # 检查输出文件
            fasta_files = [f for f in os.listdir(output_dir) if f.endswith('.fasta')]
            domain_file = os.path.join(output_dir, "domain_info.md")
            
            summary = f"""
## 执行结果
- 请求基因: {', '.join(genes)}
- 成功获取序列: {len(fasta_files)} 个
- 输出目录: {output_dir}

### 命令输出
```
{result.stdout}
```

### 运行信息
```
{result.stderr}
```
            """
            
            # 准备下载文件
            zip_path = None
            if fasta_files:
                zip_path = shutil.make_archive(
                    os.path.join(self.zip_dir, f"protein_sequences_{timestamp}"), 
                    'zip', 
                    output_dir
                )
            
            domain_content = None
            if os.path.exists(domain_file):
                with open(domain_file, 'r', encoding='utf-8') as f:
                    domain_content = f.read()
            
            return summary, zip_path, domain_content
            
        except Exception as e:
            return f"执行出错: {str(e)}", None, None
    
    def extract_sequence(self, fasta_file, start_pos: int, end_pos: int):
        """提取序列子片段"""
        if not fasta_file:
            return "请上传FASTA文件", None
        
        try:
            # 使用CommonUtils.get_output_dir生成输出目录和时间戳
            work_dir, timestamp = CommonUtils.get_output_dir("result", "extract_output", "extract")
            
            # 保存上传的文件到工作目录
            input_filename = os.path.basename(fasta_file.name)
            input_path = os.path.join(work_dir, f"input_{input_filename}")
            shutil.copy(fasta_file.name, input_path)
            
            # 执行提取命令，输出也保存在同一目录
            cmd = [
                sys.executable, "pypda.py", "seq", "extract",
                input_path, str(start_pos), str(end_pos),
                "--output_dir", work_dir
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=os.path.dirname(os.path.abspath(__file__)))
            
            # 直接将整个工作目录中的所有文件打包
            zip_path = shutil.make_archive(
                os.path.join(self.zip_dir, f"extract_results_{timestamp}"), 
                'zip', 
                work_dir
            )

            # 统计工作目录中的文件数量
            total_files = len([f for f in os.listdir(work_dir) if os.path.isfile(os.path.join(work_dir, f))])
            
            summary = f"""
## 序列提取结果
- 输入文件: {input_filename}
- 提取位置: {start_pos}-{end_pos}
- 工作目录: {work_dir}
- 总文件数: {total_files} 个

### 命令输出
```
{result.stdout}
```

### 错误信息
```
{result.stderr}
```
            """
            
            return summary, zip_path
            
        except Exception as e:
            return f"执行出错: {str(e)}", None
    
    def mutate_sequence(self, fasta_file, positions: str, amino_acids: str):
        """执行序列突变"""
        if not fasta_file:
            return "请上传FASTA文件", None
        
        try:
            # 解析输入
            pos_list = [int(p.strip()) for p in positions.split(',') if p.strip()]
            aa_list = [a.strip().upper() for a in amino_acids.split(',') if a.strip()]
            
            if len(pos_list) != len(aa_list):
                return "突变位置和氨基酸数量必须相同", None
            
            # 使用CommonUtils.get_output_dir生成输出目录和时间戳
            work_dir, timestamp = CommonUtils.get_output_dir("result", "mutate_output", "mutate")
            
            # 保存上传的文件到工作目录
            input_filename = os.path.basename(fasta_file.name)
            input_path = os.path.join(work_dir, f"input_{input_filename}")
            shutil.copy(fasta_file.name, input_path)
            
            # 执行突变命令，输出也保存在同一目录
            cmd = [
                sys.executable, "pypda.py", "seq", "mutate",
                input_path,
                "--pos"
            ] + [str(p) for p in pos_list] + ["--aa"] + aa_list + ["--output_dir", work_dir]
            
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=os.path.dirname(os.path.abspath(__file__)))
            
            # 直接将整个工作目录中的所有文件打包
            zip_path = shutil.make_archive(
                os.path.join(self.zip_dir, f"mutate_results_{timestamp}"), 
                'zip', 
                work_dir
            )

            # 统计工作目录中的文件数量
            total_files = len([f for f in os.listdir(work_dir) if os.path.isfile(os.path.join(work_dir, f))])
            
            summary = f"""
## 序列突变结果
- 输入文件: {input_filename}
- 突变位置: {positions}
- 新氨基酸: {amino_acids}
- 工作目录: {work_dir}
- 总文件数: {total_files} 个

### 命令输出
```
{result.stdout}
```

### 错误信息
```
{result.stderr}
```
            """
            
            return summary, zip_path
            
        except Exception as e:
            return f"执行出错: {str(e)}", None
    
    def align_sequences(self, fasta_files):
        """序列比对"""
        if not fasta_files or len(fasta_files) < 2:
            return "请上传至少两个FASTA文件", None
        
        try:
            # 使用CommonUtils.get_output_dir生成输出目录和时间戳
            align_temp_dir, timestamp = CommonUtils.get_output_dir("result", "align_output", "align")
            
            file_paths = []
            for i, file in enumerate(fasta_files):
                file_path = os.path.join(align_temp_dir, f"input_{i}.fasta")
                shutil.copy(file.name, file_path)
                file_paths.append(file_path)
            
            # 执行比对命令
            cmd = [sys.executable, "pypda.py", "seq", "align"] + file_paths
            
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=os.path.dirname(os.path.abspath(__file__)))
            
            # 创建zip文件包含所有输入文件和可能的输出文件
            zip_path = None
            if os.listdir(align_temp_dir):
                zip_path = shutil.make_archive(
                    os.path.join(self.zip_dir, f"align_results_{timestamp}"), 
                    'zip', 
                    align_temp_dir
                )
            
            summary = f"""
## 序列比对结果
- 输入文件数量: {len(fasta_files)}
- 临时目录: {align_temp_dir}

### 命令输出
```
{result.stdout}
```

### 错误信息
```
{result.stderr}
```
            """
            
            return summary, zip_path
            
        except Exception as e:
            return f"执行出错: {str(e)}", None
    
    def get_uniprot_data(self, protein_name: str, progress=gr.Progress()):
        """获取UniProt数据"""
        if not protein_name or not str(protein_name).strip():
            return "请输入蛋白质名称", None, None
        
        try:
            # 使用CommonUtils.get_output_dir生成输出目录和时间戳
            output_dir, timestamp = CommonUtils.get_output_dir("result", "uniprot_reports", "uniprot")
            
            # 执行命令
            cmd = [
                sys.executable, "pypda.py", "uniprot", "fetch",
                protein_name.strip(),
                "--output_dir", output_dir
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=os.path.dirname(os.path.abspath(__file__)))
            
            # 查找输出文件
            json_files = [f for f in os.listdir(output_dir) if f.endswith('.json')]
            md_files = [f for f in os.listdir(output_dir) if f.endswith('.md')]
            
            json_path = None
            md_path = None
            if json_files:
                json_path = os.path.join(output_dir, json_files[0])
            if md_files:
                md_path = os.path.join(output_dir, md_files[0])
            
            # 创建zip文件
            zip_path = None
            if json_files or md_files:
                zip_path = shutil.make_archive(
                    os.path.join(self.zip_dir, f"uniprot_results_{timestamp}"), 
                    'zip', 
                    output_dir
                )
            
            summary = f"""
## UniProt数据获取结果
- 蛋白质名称: {protein_name}
- 输出目录: {output_dir}
- JSON文件: {len(json_files)} 个
- Markdown文件: {len(md_files)} 个

### 命令输出
```
{result.stdout}
```

### 错误信息
```
{result.stderr}
```
            """
            
            # 读取markdown内容
            md_content = None
            if md_path and os.path.exists(md_path):
                with open(md_path, 'r', encoding='utf-8') as f:
                    md_content = f.read()
            
            return summary, zip_path, md_content
            
        except Exception as e:
            return f"执行出错: {str(e)}", None, None
    
    def process_pdb(self, protein_name: str):
        """处理PDB文件"""
        if not protein_name or not str(protein_name).strip():
            return "请输入蛋白质名称", None
        
        try:
            # 使用CommonUtils.get_output_dir生成输出目录和时间戳
            output_dir, timestamp = CommonUtils.get_output_dir("result", "pdb_output", "pdb")
            
            # 执行命令
            cmd = [
                sys.executable, "pypda.py", "pdb",
                protein_name.strip(),
                output_dir
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=os.path.dirname(os.path.abspath(__file__)))
            
            # 创建zip文件
            zip_path = None
            if os.listdir(output_dir):
                zip_path = shutil.make_archive(
                    os.path.join(self.zip_dir, f"pdb_results_{timestamp}"), 
                    'zip', 
                    output_dir
                )
            
            summary = f"""
## PDB处理结果
- 蛋白质名称: {protein_name}
- 输出目录: {output_dir}
- 文件数量: {len(os.listdir(output_dir))}

### 命令输出
```
{result.stdout}
```

### 错误信息
```
{result.stderr}
```
            """
            
            return summary, zip_path
            
        except Exception as e:
            return f"执行出错: {str(e)}", None

# 创建Gradio界面
def create_interface():
    """创建Gradio界面"""
    interface = PyPDAGradioInterface()
    
    with gr.Blocks(title="PyPDA - 蛋白质数据分析工具", theme=gr.themes.Soft()) as app:
        gr.Markdown("""
        # 🧬 PyPDA - 蛋白质数据分析工具
        
        欢迎使用PyPDA的Web界面！这个工具提供了多种蛋白质数据分析功能，包括：
        
        - 🎯 基因序列获取和结构域分析
        - ✂️ 序列提取和突变
        - 🔄 序列比对
        - 📊 UniProt数据获取和分析
        - 🏗️ PDB文件处理
        """)
        
        with gr.Tab("🎯 获取蛋白质序列"):
            with gr.Row():
                with gr.Column():
                    gr.Markdown("### 获取蛋白质序列和结构域信息\n输入基因名称列表，获取对应的蛋白质序列和结构域信息\n")
                    genes_input = gr.Textbox(
                        label="基因名称",
                        placeholder="输入基因名称，用空格分隔 (例如: BRCA1 TP53 EGFR)",
                        lines=3
                    )
                    get_sequences_btn = gr.Button("获取序列", variant="primary")
                
                with gr.Column():
                    sequences_summary = gr.Markdown()
                    sequences_download = gr.File(label="下载序列文件")
                    domain_info = gr.Markdown()
        
        with gr.Tab("✂️ 序列提取"):
            with gr.Row():
                with gr.Column():
                    gr.Markdown("### 提取序列子片段\n从FASTA文件中提取指定位置的子序列\n")
                    extract_file = gr.File(
                        label="上传FASTA文件",
                        file_types=[".fasta", ".fa", ".txt"]
                    )
                    extract_start = gr.Number(label="起始位置", value=1, minimum=1)
                    extract_end = gr.Number(label="结束位置", value=100, minimum=1)
                    extract_btn = gr.Button("提取序列", variant="primary")
                
                with gr.Column():
                    extract_summary = gr.Markdown()
                    extract_download = gr.File(label="下载提取结果")
        
        with gr.Tab("🧬 序列突变"):
            with gr.Row():
                with gr.Column():
                    gr.Markdown("### 执行序列突变\n对FASTA文件中的序列进行指定位置的点突变\n")
                    mutate_file = gr.File(
                        label="上传FASTA文件",
                        file_types=[".fasta", ".fa", ".txt"]
                    )
                    mutate_positions = gr.Textbox(
                        label="突变位置",
                        placeholder="用逗号分隔的位置 (例如: 10,20,30)"
                    )
                    mutate_amino_acids = gr.Textbox(
                        label="新氨基酸",
                        placeholder="用逗号分隔的氨基酸 (例如: A,G,L)"
                    )
                    mutate_btn = gr.Button("执行突变", variant="primary")
                
                with gr.Column():
                    mutate_summary = gr.Markdown()
                    mutate_download = gr.File(label="下载突变结果")
        
        with gr.Tab("🔄 序列比对"):
            with gr.Row():
                with gr.Column():
                    gr.Markdown("### 序列比对\n对多个FASTA文件进行序列比对\n")
                    align_files = gr.File(
                        label="上传FASTA文件",
                        file_types=[".fasta", ".fa", ".txt"],
                        file_count="multiple"
                    )
                    align_btn = gr.Button("开始比对", variant="primary")
                
                with gr.Column():
                    align_summary = gr.Markdown()
                    align_download = gr.File(label="下载比对结果")
        
        with gr.Tab("📊 UniProt数据"):
            with gr.Row():
                with gr.Column():
                    gr.Markdown("### 获取UniProt数据并生成报告\n输入基因名称，获取详细的UniProt数据报告\n")
                    uniprot_input = gr.Textbox(
                        label="蛋白质名称",
                        placeholder="输入蛋白质名称 (例如: TP53, BRCA1)"
                    )
                    get_uniprot_btn = gr.Button("获取数据", variant="primary")
                
                with gr.Column():
                    uniprot_summary = gr.Markdown()
                    uniprot_zip = gr.File(label="下载结果文件")
                    uniprot_report = gr.Markdown()
        
        with gr.Tab("🏗️ PDB处理"):
            with gr.Row():
                with gr.Column():
                    gr.Markdown("### PDB文件处理\n下载和处理PDB结构文件\n")
                    pdb_protein_input = gr.Textbox(
                        label="蛋白质名称",
                        placeholder="输入蛋白质名称 (例如: BRCA1, TP53)"
                    )
                    process_pdb_btn = gr.Button("处理PDB", variant="primary")
                
                with gr.Column():
                    pdb_summary = gr.Markdown()
                    pdb_download = gr.File(label="下载PDB结果")
        
        gr.Markdown("""
        ---
        ### 📋 使用说明

        1. **获取蛋白质序列**: 输入基因名称（如 BRCA1 TP53），获取对应的蛋白质序列和结构域信息
        2. **序列提取**: 上传FASTA文件，指定起始和结束位置提取子序列
        3. **序列突变**: 上传FASTA文件，指定突变位置和新的氨基酸
        4. **序列比对**: 上传多个FASTA文件进行两两比对
        5. **UniProt数据**: 获取指定蛋白质的详细UniProt信息
        6. **PDB处理**: 下载和处理PDB结构文件

        ### ⚠️ 注意事项

        - 所有操作结果都会保存在 result 目录中，请及时下载需要的文件
        - 处理大文件可能需要一些时间，请耐心等待
        - 如有错误，请查看命令输出信息
        """)
        
        # 绑定事件
        get_sequences_btn.click(
            interface.get_protein_sequences,
            inputs=[genes_input],
            outputs=[sequences_summary, sequences_download, domain_info]
        )
        
        extract_btn.click(
            interface.extract_sequence,
            inputs=[extract_file, extract_start, extract_end],
            outputs=[extract_summary, extract_download]
        )
        
        mutate_btn.click(
            interface.mutate_sequence,
            inputs=[mutate_file, mutate_positions, mutate_amino_acids],
            outputs=[mutate_summary, mutate_download]
        )
        
        align_btn.click(
            interface.align_sequences,
            inputs=[align_files],
            outputs=[align_summary, align_download]
        )
        
        get_uniprot_btn.click(
            interface.get_uniprot_data,
            inputs=[uniprot_input],
            outputs=[uniprot_summary, uniprot_zip, uniprot_report]
        )
        
        process_pdb_btn.click(
            interface.process_pdb,
            inputs=[pdb_protein_input],
            outputs=[pdb_summary, pdb_download]
        )
    
    return app

if __name__ == "__main__":
    import argparse
    
    # 解析命令行参数
    parser = argparse.ArgumentParser(description='PyPDA Gradio Web界面')
    parser.add_argument('--port', type=int, default=7860, help='服务器端口，默认7860')
    args = parser.parse_args()
    port = args.port
    
    app = create_interface()
    app.launch(
        server_name="0.0.0.0",
        server_port=port,
        share=False,
        debug=True
    )