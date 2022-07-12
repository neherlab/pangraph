from setuptools import find_packages, setup

setup(
    name="pypangraph",
    version="0.1",
    description="package to manipulate pangenomes produced by PanGraph",
    author="Marco Molari",
    author_email="mmolari91@gmail.com",
    url="https://github.com/mmolari/pypangraph",
    packages=find_packages(),
    install_requires = [
        "matplotlib",
        "numpy",
        "pandas",
        "biopython",
        "scipy",
        "networkx",
    ]
)
