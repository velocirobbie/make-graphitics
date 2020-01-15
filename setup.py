from setuptools import setup, find_packages

setup(
    name="makegraphitics",
    version="0.2",
    description="""Library to build graphene and graphite based structures
                   for atomistic simulation""",
    url="https://github.com/velocirobbie/make-graphitics",
    author="Robert C Sinclair",
    packages=find_packages(),
    include_package_data=True
)
