from setuptools import setup, find_packages

setup(
    name='automatchnames',
    version='1.3.1',
    packages=find_packages(),
    package_data={"wcvp_download": ["inputs/*", "inputs/wgsrpd-master/level3/*"]},
    install_requires=[
        'pandas>=2.1.4',
        'numpy>=1.26',
        'requests>=2.31',
        'tqdm>=4.66',
        'typing',
        'zstandard'
    ],
    extras_require={
        'dist_plots': ["matplotlib", 'cartopy', 'fiona', 'pillow']
    },
    url='https://github.com/alrichardbollans/automatchnames',
    license='GNU v.3',
    author='Adam Richard-Bollans',
    description='A package for downloading the WCVP and matching names to it',
    long_description=open('readme.md', encoding="utf8").read()
)
