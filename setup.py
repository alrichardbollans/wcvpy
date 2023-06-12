from setuptools import setup, find_packages

setup(
    name='automatchnames',
    version='1.2.0',
    packages=find_packages(),
    install_requires=[
        'pandas==2.0.1',
        'numpy==1.24.2',
        'requests==2.28.2',
        'tqdm==4.64.1',
        'typing',
        'hashlib'

    ],
    url='https://github.com/alrichardbollans/automatchnames',
    license='GNU v.3',
    author='Adam Richard-Bollans',
    description='A package for downloading WCVP and matching names to it',
    long_description=open('readme.md', encoding="utf8").read()
)
