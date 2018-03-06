from setuptools import setup, find_packages

setup(
    name="gocamgen",
    version='0.0.1',
    # packages=["gocamgen"],
    packages=find_packages(),
    author="dustine32",
    author_email="debert@usc.edu",
    description="Python library for constructing GO-CAM model RDF",
    long_description=open("README.md").read(),
    url="https://github.com/dustine32/gocamgen",
    install_requires=[
        "ontobio"
    ]
)