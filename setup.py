from setuptools import setup, find_packages

setup(
    name='smog3',
    version='2.6-beta',
    packages=find_packages(where='.'),
    package_dir={'': '.'},
    include_package_data=True,
    package_data={
        'smog3': ['share/templates/**/*'],
    },
    install_requires=[
        'numpy',
    ],
    entry_points={
        'console_scripts': [
            'smog3=smog3.src.smog3:main',
        ],
    },
    author='SMOG Team',
    description='Structure-based Model (SMOG) software in Python',
    url='http://smog-server.org',
)
