import glob
from setuptools import setup, find_packages

try:
    import multiprocessing
except ImportError:
    pass

setup(
    name='eupathws',
    version='0.1',
    description='Python interface for exporting EuPathDB data from web services',
    packages = find_packages(),
    author='Sascha Steinbiss',
    author_email='sascha@steinbiss.name',
    url='https://github.com/satta/eupathws',
    scripts=glob.glob('scripts/*'),
    license='ISC',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
)
