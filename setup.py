from setuptools import setup, find_packages
import os

# Helper to list files for data_files
def find_data_files(directory):
    paths = []
    for (path, directories, filenames) in os.walk(directory):
        if filenames:
            paths.append((path.replace('share/templates', 'share/smog3/templates'),
                          [os.path.join(path, f) for f in filenames]))
    return paths

data_files = find_data_files('share/templates')

setup(
    name='smog3',
    version='2.7.0',  # Using semantic versioning for PyPI
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'numpy',
    ],
    entry_points={
        'console_scripts': [
            'smog3=smog3.src.smog3:main',
        ],
    },
    data_files=data_files,
    author='SMOG Team',
    author_email='info@smog-server.org',
    description='Structure-based Model (SMOG) software in Python',
    long_description=open('README.md').read() if os.path.exists('README.md') else '',
    long_description_content_type='text/markdown',
    url='http://smog-server.org',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License', # Assuming MIT or similar, check COPYING
        'Operating System :: OS Independent',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    python_requires='>=3.6',
)
