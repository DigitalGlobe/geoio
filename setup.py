from setuptools import setup, find_packages
import os

# import __version__
exec(open('geoio/_version.py').read())

install_requires = [
    'gdal',
    'xmltodict',
    'pytz',
    'tzwhere',
    'ephem',
    'numpy',
    'tinytools'
    ]
    # optional libraries for additional capabilities:
    # 'matplotlib'
    # 'cv2'
    # 'numba'

tests_require = [
    'dgsamples'
]

setup(
    name='geoio',
    version=__version__,
    author='Nathan Longbotham',
    author_email='nlongbotham@digitalglobe.com',
    packages=find_packages(),
    description='Geo image reading/writing tools.',
    long_description=open('README.rst').read(),
    install_requires=install_requires,
    tests_require = tests_require,
    scripts=['bin/run_spectral_conversions.py','bin/find_files.py']
)
