import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()
    
setuptools.setup(
    name="apogeesb2",
    version="1.2",
    author="Marina Kounkel",
    author_email="marina.kounkel@wwu.edu",
    description="Examine CCFs in the APOGEE apstar files to search for double lined spectroscopic binaries",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/mkounkel/apogeesb2",
    packages=['apogeesb2'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=['gausspy','astropy',],
    entry_points={'console_scripts': ['apogeesb2=apogeesb2.apogeesb2:run']},
    python_requires='>=3.6',
)
