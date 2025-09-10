from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="conflux",
    version="1.1.1",
    author="Xianyi Zhang",
    author_email="zhang39@llnl.gov",
    description="A package to calculate neutrino flux from beta decaying sources",
    long_description="README.md",
    long_description_content_type="text/markdown",
    url="https://github.com/CNFLUX/conflux",
    project_urls={
        "Bug Tracker": "https://github.com/CNFLUX/conflux",
    },

    packages=find_packages(),
    package_data={
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
	install_requires=['numpy', 'scipy>=1.8.1', 'tqdm', 'matplotlib', 'iminuit', 'fortranformat', 'pandas', 'xraydb'],
    python_requires=">=3.6",
)
