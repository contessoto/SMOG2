from setuptools import setup, find_packages

setup(
    name='smog3',
    version='2.7beta',
    packages=find_packages(where='.'),
    package_dir={'': '.'},
    include_package_data=True,
    package_data={
        # This only works if files are inside the package. share is outside.
    },
    data_files=[
        ('share/smog3/templates', [
            'share/templates/ff.info'
        ]),
        ('share/smog3/templates/SBM_AA', [
            'share/templates/SBM_AA/AA-whitford09.b',
            'share/templates/SBM_AA/AA-whitford09.sif',
            'share/templates/SBM_AA/default.map',
            'share/templates/SBM_AA/.citation',
            'share/templates/SBM_AA/AA-whitford09.bif',
            'share/templates/SBM_AA/AA-whitford09.nb'
        ]),
        ('share/smog3/templates/SBM_AA+gaussian', [
            'share/templates/SBM_AA+gaussian/AA+gaussian-noel12.b',
            'share/templates/SBM_AA+gaussian/AA+gaussian-noel12.sif',
            'share/templates/SBM_AA+gaussian/AA+gaussian-noel12.nb',
            'share/templates/SBM_AA+gaussian/default.map',
            'share/templates/SBM_AA+gaussian/AA+gaussian-noel12.bif',
            'share/templates/SBM_AA+gaussian/.citation'
        ]),
        ('share/smog3/templates/SBM_calpha', [
            'share/templates/SBM_calpha/CA-clementi00.nb',
            'share/templates/SBM_calpha/default.map',
            'share/templates/SBM_calpha/CA-clementi00.sif',
            'share/templates/SBM_calpha/.citation',
            'share/templates/SBM_calpha/CA-clementi00.b',
            'share/templates/SBM_calpha/CA-clementi00.bif'
        ]),
        ('share/smog3/templates/SBM_calpha+gaussian', [
            'share/templates/SBM_calpha+gaussian/CA+gaussian-lammert09.nb',
            'share/templates/SBM_calpha+gaussian/default.map',
            'share/templates/SBM_calpha+gaussian/.citation',
            'share/templates/SBM_calpha+gaussian/CA+gaussian-lammert09.b',
            'share/templates/SBM_calpha+gaussian/CA+gaussian-lammert09.sif',
            'share/templates/SBM_calpha+gaussian/CA+gaussian-lammert09.bif'
        ]),
    ],
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
