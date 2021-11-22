import pathlib
from setuptools import setup

# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / "README.md").read_text()

# This call to setup() does all the work
setup(
    name="phylotreelib",
    version="1.4.0",
    description="Analyze and manipulate phylogenetic trees",
    long_description=README,
    long_description_content_type="text/markdown",
    url="https://github.com/agormp/phylotreelib",
    author="Anders Gorm Pedersen",
    author_email="agpe@dtu.dk",
    license="GPLv3+",
    classifiers=[
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
    ],
    py_modules=['phylotreelib'],
    install_requires=["numpy"],
)