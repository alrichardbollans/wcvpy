from setuptools import setup

setup(
    name='automatchnames',
    version='0.0.2',
    packages=['wcvp_name_matching', 'wcvp_download'],
    install_requires=[
        'pandas==1.4.1',
        'numpy~=1.22.1',
        'requests~=2.27.1',
        'tqdm~=4.62.3',
        'typing~=3.7.4.3'

    ],
    url='https://github.com/alrichardbollans/automatchnames',
    license='GNU v.3',
    author='Adam Richard-Bollans',
    description='A package for automating name matching of scientific names',
    long_description=open('readme.md', encoding="utf8").read()
)
