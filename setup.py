from setuptools import setup, find_packages
import os

# Set version number
vpath = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(vpath,'VERSION.txt')) as f:
    VERSION = f.read().strip('\n')

install_requires = [
    'gdal',
    'xmltodict',
    'pytz',
    'tzwhere',
    'ephem',
    'numpy',
    'tinytools'
    ]

tests_require = [
    'dgsamples'
]

setup(
    name='geoio',
    version=VERSION,
    author='Nathan Longbotham',
    author_email='nlongbotham@digitalglobe.com',
    packages=find_packages(),
    description='Geo image reading/writing tools.',
    long_description=open('README.rst').read(),
    install_requires=install_requires,
    tests_require = tests_require,
    scripts=['bin/run_spectral_conversions.py','bin/find_files.py']
)
