from setuptools import setup, find_packages

with open('README.md', 'r') as fh:
    long_description = fh.read()


setup(
    name='gendock',
    version='1.0',
    description='Python package for generating and docking molecules against macromolecules',
    url='https://github.com/philspence/gendock',
    author='Phil Spence',
    author_email='philspence91@gmail.com',
    packages=find_packages(),
    long_description=long_description,
    long_description_content_type='text/markdown',
    python_requires='>=3.7',
    classifiers=[
        'Programming Language :: Python :: 3',
        'Operating System :: MacOS',
        'Operating System :: POSIX :: Linux'
    ]
)
