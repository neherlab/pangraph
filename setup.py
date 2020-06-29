from pathlib import Path
import setuptools
import sys

min_version = (3, 6)

if sys.version_info < min_version:
    error = """
Python {0} or above is required.
This may be due to an out of date pip.
Make sure you have pip >= 9.0.1.
""".format('.'.join(str(n) for n in min_version)),
    sys.exit(error)

base_dir = Path(__file__).parent.resolve()
version_file = base_dir / "pangraph/__version__.py"
readme_file = base_dir / "README.md"

# Eval the version file to get __version__; avoids importing our own package
with version_file.open() as f:
    exec(f.read())

setuptools.setup(
    name = "pangraph",
    version = __version__,
    author = "Pangraph developers",
    author_email = "nbnoll@eml.cc, richard.neher@unibas.ch",
    description = "A bioinformatics toolkit for pangenomic analysis",
    keywords = "pangraph, molecular epidemiology",
    url = "https://github.com/nnoll/pangraph",
    project_urls = {
        "Bug Reports": "https://github.com/nnoll/pangraph/issues",
        "Source": "https://github.com/nnoll/pangraph",
    },
    packages = setuptools.find_packages(),
    package_data = {'pangraph': ['data/*']},
    data_files = [("", ["LICENSE"])],
    python_requires = f">={min_version[0]}.{min_version[1]}",
    install_requires = [
        "numpy",
        "biopython",
        "portion",
        "scikit-bio",
        "cigar",
        "suffix-tree",
    ],
    extras_require = {
        'dev': [
            "ipython",
            "matplotlib",
            "ipdb",
         ],
    },
    classifiers = [
        "Topic :: Scientific/Engineering :: Bio-Informatics",

        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
    ],
    entry_points = {
        "console_scripts": [
            "pangraph = pangraph.__main__:main",
        ]
    }
)
