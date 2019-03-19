import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="ssdb",
    version="0.0.1",
    author="Robin Ekelund, Simon Pfreundschuh",
    author_email="simon.pfreundschuh@chalmers.se",
    description="Python interface for the ARTS single scattering database.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/simonpf/ssdb",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
