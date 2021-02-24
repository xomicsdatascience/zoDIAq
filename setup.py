import setuptools


def readfile(filename):
    with open(filename, 'r+') as f:
        return f.read()

setuptools.setup(
    name="csodiaq", # Replace with your own username
    version="0.0.1",
    author="Caleb Cranney",
    author_email="caleb.cranney.app@gmail.com",
    description="CsoDIAq Package",
    long_description=readfile('README.ipynb'),
    long_description_content_type="text/markdown",
    url="https://github.com/CCranney/csoDIAq",
    classifiers=[
        'Intended Audience :: Science/Research',
        "Programming Language :: Python :: 3",
    ],
    py_modules=['csodiaq','csodiaq_menu_functions','csodiaq_base_functions','csodiaq_gui','idpicker'],
    packages=setuptools.find_packages(),
    python_requires='>=3.6',
    license=readfile('LICENSE'),
    entry_points={
        'console_scripts': [
            'csodiaq = csodiaq:main'
        ]
        #gui_scripts?
    },
    #install_requires=[
    #'certifi==2020.12.5',
    #'cycler==0.10.0',
    #'kiwisolver==1.3.1',
    #'lxml==4.6.2',
    #'matplotlib==3.3.4',
    #'numpy==1.20.1',
    #'pandas==1.2.2',
    #'Pillow==8.1.0',
    #'pyparsing==2.4.7',
    #'PyQt5==5.15.2',
    #'PyQt5-sip==12.8.1',
    #'pyteomics==4.4.1',
    #'python-dateutil==2.8.1',
    #'pytz==2021.1',
    #'six==1.15.0',
    #],
)
