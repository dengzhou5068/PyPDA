#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PyPDA Gradio启动脚本
一键启动PyPDA的Web界面
"""

import os
import sys
import subprocess
import webbrowser
import time
from pathlib import Path
import argparse

def check_dependencies():
    """检查必要的依赖"""
    try:
        import gradio
        print("✅ Gradio 已安装")
    except ImportError:
        print("❌ Gradio 未安装，正在安装...")
        subprocess.run([sys.executable, "-m", "pip", "install", "gradio>=4.0.0"])
        print("✅ Gradio 安装完成")

def main():
    """主函数"""
    # 解析命令行参数
    parser = argparse.ArgumentParser(description='启动PyPDA Web界面')
    parser.add_argument('--port', type=int, default=7860, help='服务器端口，默认7860')
    args = parser.parse_args()
    port = args.port
    
    print("🧬 启动 PyPDA Web界面...")
    print("=" * 50)
    
    # 检查依赖
    check_dependencies()
    
    # 获取当前目录
    current_dir = Path(__file__).parent.absolute()
    os.chdir(current_dir)
    
    print(f"📁 工作目录: {current_dir}")
    print(f"🚀 正在启动Gradio应用，端口: {port}...")
    
    try:
        # 启动Gradio应用
        from gradio_app import create_interface
        app = create_interface()
        
        print("\n🌐 Web界面启动成功!")
        print("=" * 50)
        print(f"访问地址: http://localhost:{port}")
        print("\n📋 功能说明:")
        print("  - 🎯 获取蛋白质序列和结构域信息")
        print("  - ✂️ 序列提取和突变")
        print("  - 🔄 序列比对")
        print("  - 📊 UniProt数据获取和分析")
        print("  - 🏗️ PDB文件处理")
        print("\n📁 文件存储结构:")
        print("   📂 result/")
        print("      📂 protein_sequences/ - 蛋白质序列结果")
        print("      📂 extract_output/ - 序列提取结果")
        print("      📂 mutate_output/ - 序列突变结果")
        print("      📂 align_output/ - 序列比对结果")
        print("      📂 uniprot_reports/ - UniProt数据报告")
        print("      📂 pdb_output/ - PDB文件处理结果")
        print("   📌 每个功能按时间戳创建子文件夹，避免文件混淆")
        print("\n⚠️  按 Ctrl+C 停止服务")
        print("=" * 50)
        
        # 自动打开浏览器
        try:
            time.sleep(2)
            webbrowser.open(f"http://localhost:{port}")
        except:
            pass
        
        # 启动应用
        app.launch(
            server_name="0.0.0.0",
            server_port=port,
            share=False,
            debug=False,
            show_error=True
        )
        
    except KeyboardInterrupt:
        print("\n👋 正在关闭服务...")
    except Exception as e:
        print(f"❌ 启动失败: {e}")
        print("请检查依赖是否安装正确")

if __name__ == "__main__":
    main()