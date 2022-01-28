from setuptools import setup, Extension
from os import path
import io


this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, "README.md"), encoding='utf-8') as f:
    long_description = f.read()


if __name__ == "__main__":
    setup(
        name="melt",
        long_description=long_description,
        long_description_content_type="text/markdown",
        author="Ben Mather",
        author_email="brmather1@gmail.com",
        url="https://github.com/brmather/melt",
        version='0.1',
        description="Python tool for computing wet melting",
        install_requires=["numpy", "scipy>=0.11.0"],
        python_requires=">=2.7, >=3.5",
        packages=["melt"],
        package_data={
            "melt": [
                "Notebooks/*.ipynb",
            ]
        },
        classifiers=[
            "Programming Language :: Python :: 2",
            "Programming Language :: Python :: 2.7",
            "Programming Language :: Python :: 3",
            "Programming Language :: Python :: 3.5",
            "Programming Language :: Python :: 3.6",
            "Programming Language :: Python :: 3.7",
            "Programming Language :: Python :: 3.8",
            "Programming Language :: Python :: 3.9",
            "Programming Language :: Python :: 3.10",
        ],
    )