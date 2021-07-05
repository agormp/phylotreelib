import pathlib
from setuptools import setup

# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / "README.md").read_text()

# This call to setup() does all the work
setup(
    name="phylotreelib",
    version="1.3.1",
    description="Analyze and manipulate phylogenetic trees",
    long_description=README,
    long_description_content_type="text/markdown",
    url="https://github.com/agormp/phylotreelib",
    author="Anders Gorm Pedersen",
    author_email="agpe@dtu.dk",
    license="MIT",
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
    ],
    py_modules=['phylotreelib'],
    install_requires=["numpy"],
)