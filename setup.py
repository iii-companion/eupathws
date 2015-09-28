import glob
from setuptools import setup, find_packages

setup(
    name='eupathtables',
    version='0.1',
    description='Python interface for reading and converting EuPathDB flat file dumps',
    packages = find_packages(),
    author='Sascha Steinbiss',
    author_email='ss34@sanger.ac.uk',
    url='https://github.com/satta/eupathtables',
    scripts=glob.glob('scripts/*'),
    license='GPLv3',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
)
