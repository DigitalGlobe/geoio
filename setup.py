from setuptools import setup, find_packages

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
    version='1.1.10',
    author='Nathan Longbotham',
    author_email='nlongbotham@digitalglobe.com',
    packages=find_packages(),
    description='Geo image reading/writing tools.',
    long_description=open('README.rst').read(),
    install_requires=install_requires,
    tests_require = tests_require,
    scripts=['bin/run_spectral_conversions.py','bin/find_files.py']
)
