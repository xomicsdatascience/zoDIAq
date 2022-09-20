import setuptools
import os

def readfile(filename):
    with open(filename, 'r+') as f:
        return f.read()


def get_version(module_path: str) -> str:
    """Gets the version from the specified file. Expected to be of format __version__ = 'x.y.z'"""
    f = open(module_path, 'r')
    lines = f.read().splitlines()
    for l in lines:
        if l.startswith('__version__'):
            # Found version; parse
            return l.split('=')[-1].strip().strip("'").strip('"')

pwd = os.path.abspath(os.path.dirname(__file__))
version = get_version(os.path.join(pwd, 'csodiaq/__init__.py'))

setuptools.setup(
    name="csodiaq",
    version=version,
    author="Meyer Lab",
    author_email="",
    description="CsoDIAq Package",
    long_description=readfile('README.md'),
    long_description_content_type="text/markdown",
    url="",
    classifiers=[
        "Development Status :: 4 - Beta",
        'Intended Audience :: Science/Research',
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License"
    ],
    packages=setuptools.find_packages(),
    python_requires='>=3.8',
    license="MIT",
    entry_points={
        'console_scripts': [
            'csodiaq = csodiaq.csodiaq:main'
        ]
    },
    install_requires=[
        'pyteomics>=4.4.1',
        'matplotlib>=3.3.4',
        'numba>=0.53.1',
        'numpy>=1.20.1',
        'pandas>=1.2.2',
        'Bio>=0.4.1',
        'PyQt5>=5.15.4',
        'lxml>=4.6.2',
        'scikit-learn>=1.1.2'
    ],
)
