from setuptools import setup, find_packages
import os

# 读取README.md作为长描述
with open(os.path.join(os.path.dirname(__file__), 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

# 读取requirements.txt作为依赖列表
with open(os.path.join(os.path.dirname(__file__), 'requirements.txt'), encoding='utf-8') as f:
    requirements = [line.strip() for line in f.readlines() if line.strip() and not line.startswith('#')]

setup(
    name='pypda',
    version='0.2.0',
    description='Protein Domain Analysis Toolkit - 蛋白质结构域分析工具包',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Zhou Deng',
    author_email='dengzho5068@foxmail.com',
    url='https://gitee.com/coding_playground/pypda',
    packages=find_packages(),
    package_data={
        '': ['exclude_residues.ini'],
    },
    include_package_data=True,
    install_requires=requirements,
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.8',
    keywords=['bioinformatics', 'protein analysis', 'domain prediction', 'sequence analysis'],
    project_urls={
        'Bug Reports': 'https://gitee.com/coding_playground/py-pda/issues',
        'Source': 'https://gitee.com/coding_playground/py-pda',
    },
)