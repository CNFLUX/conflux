from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="conflux",
    version="1.1.3",
    author="Xianyi Zhang",
    author_email="zhang39@llnl.gov",
    description="A package to calculate neutrino flux from beta decaying sources",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/CNFLUX/conflux",
    project_urls={
        "Bug Tracker": "https://github.com/CNFLUX/conflux",
    },

    packages=find_packages(),
    # Use MANIFEST.in to control what data files are included
    # include_package_data=True means: include everything in MANIFEST.in
    include_package_data=True,
    # Explicitly exclude covariance matrices (too large - 3.5 GB)
    exclude_package_data={
        '': ['*cov*.csv', '*corr*.csv'],
    },
    entry_points={
        'console_scripts': [
            'conflux-setup=conflux.cli_setup:main',
            'conflux-update-endf=update_endf_database:main',
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Physics",
        "License :: OSI Approved :: MIT License",
    ],
    install_requires=['numpy', 'scipy>=1.8.1', 'tqdm', 'matplotlib', 'iminuit', 'fortranformat', 'pandas', 'xraydb'],
    python_requires=">=3.6",
)
