#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PyPDA 命令行入口脚本
"""

import sys
import os

# 将当前目录添加到Python路径
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from main import PypdaApp


if __name__ == "__main__":
    app = PypdaApp()
    app.run()
