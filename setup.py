from io import open

from setuptools import find_packages, setup

# following src dir layout according to
# https://blog.ionelmc.ro/2014/05/25/python-packaging/#the-structure
version = "0.0.10"
setup(
    name="varianteval",
    version=version,
    description="Library to help compare VCF files",
    url="https://broadinstitute.github.io/varianteval/",
    author="Kiran V Garimella",
    author_email="kiran@broadinstitute.org",
    license="BSD 3-Clause",
    long_description=open("README.rst").read(),
    install_requires=[
        'cython',
        'mypy',
        'click',
        'click-log',
        'pandas',
        'numpy',
        'matplotlib',
        'pysam',
        'tqdm',
    ],
    tests_require=["coverage", "pytest"],
    python_requires=">=3.6",
    packages=find_packages("varianteval"),
    package_dir={"": "varianteval"},
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Natural Language :: English",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: Implementation :: CPython",
    ],
    entry_points={"console_scripts": ["varianteval=varianteval.__main__:main_entry"]},
    include_package_data=True,
)
