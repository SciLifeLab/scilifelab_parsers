from setuptools import setup, find_packages

try:
    with open("requirements.txt", "r") as f:
        install_requires = [x.strip() for x in f.readlines()]
except IOError:
    install_requires = []

setup(name = "scilifelab_parsers",
    version = "1.0",
    author = "Maya Brandi",
    author_email = "maya.brandi@scilifelab.se",
    description = "Package for scilife parser modules",
    install_requires = install_requires,
    packages = find_packages()
    )
