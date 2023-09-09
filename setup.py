from setuptools import find_packages, setup

with open("Readme.md", mode="r", encoding="utf-8") as f:
    long_description = f.read()

setup(
    name="ABInitioREMC",
    version="0.0.10",
    description="An implementation of the Replica Exchange Monte Carlo (REMC) method to resolve the AB Initio problem (NP-Hard) using Python.",
    package_dir={"": "app"},
    packages=find_packages(where="app"),
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/lounesmechouek/ABInitioREMC",
    author="Lounes Mechouek",
    author_email="contact@mechouek.com",
    license="MIT",
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.10.11",
        "Operating System :: OS Independent",
    ],
    install_requires=["numpy >= 1.25.2", "plotly >= 5.16.1", "streamlit >= 1.26.0"],
    extras_require={
        "dev": [
            "mkdocs >= 1.5.2",
            "mkdocs-material >= 9.2.8",
            "mkdocstrings >= 0.23.0",
            "mkdocstrings-python >= 1.6.2",
        ],
    },
    python_requires=">=3.10.11",
)
