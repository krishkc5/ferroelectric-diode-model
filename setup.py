from setuptools import setup, find_packages

setup(
    name="fed_model",
    version="0.1.0",
    description="Ferroelectric Diode Physics Model",
    packages=find_packages(),
    install_requires=[
        "numpy>=1.21.0",
        "scipy>=1.7.0", 
        "matplotlib>=3.5.0",
        "pandas>=1.3.0",
    ],
    python_requires=">=3.8",
)
