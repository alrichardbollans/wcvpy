from setuptools import setup, find_packages

setup(
    name='wcvpy',
    version='1.3.3',
    packages=find_packages(),
    package_data={"wcvpy": ["wcvp_download/inputs/*", "wcvp_download/inputs/wgsrpd-master/level3/*"]},
    install_requires=[
        'pandas',
        'numpy',
        'requests',
        'tqdm',
        'typing',
        'zstandard'
    ],
    extras_require={
        'dist_plots': ["matplotlib", 'cartopy', 'fiona', 'pillow']
    },
    url='https://github.com/alrichardbollans/wcvpy',
    license='GNU v.3',
    author='Adam Richard-Bollans',
    description='A package for downloading the WCVP and matching names to it',
    long_description=open('readme.md', encoding="utf8").read()
)
