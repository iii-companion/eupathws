import glob
from setuptools import setup, find_packages

try:
    import multiprocessing
except ImportError:
    pass

setup(
    name='eupathtables',
    version='0.1',
    description='Python interface for exporting EuPathDB data from web services',
    packages = find_packages(),
    author='Sascha Steinbiss',
    author_email='ss34@sanger.ac.uk',
    url='https://github.com/satta/eupathtables',
    test_suite='nose.collector',
    tests_require=['nose >= 1.3'],
    scripts=glob.glob('scripts/*'),
    license='ISC',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
)
