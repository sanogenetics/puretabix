from setuptools import setup

setup(
    name="puretabix",
    version="5.4.0",
    author="Adam Faulconbridge",
    author_email="afaulconbridge@googlemail.com",
    packages=["puretabix"],
    package_data={"puretabix": ["py.typed"]},
    description="Pure python implementation Tabix reader.",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/sanogenetics/puretabix",
    install_requires=[],
    extras_require={
        "dev": [
            "pytest-cov",
            "flake8",
            "black",
            "mypy",
            "pylint",
            "pip-tools",
            "pipdeptree",
            "pre-commit",
            "twine",
            "scalene",
        ],
    },
)
