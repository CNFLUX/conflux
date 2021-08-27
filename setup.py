#FIXME: edit the setup to match this package

import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="conflux",
    version="0.0.1",
    author="Xianyi Zhang",
    author_email="zhang39@llnl.gov",
    description="A package to calculate neutrino flux from beta decaying sources",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://lc.llnl.gov/bitbucket/projects/CFX/repos/conflux/browse",
    project_urls={
        "Bug Tracker": "https://lc.llnl.gov/bitbucket/projects/CFX/repos/conflux/browse",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    python_requires=">=3.6",
)
