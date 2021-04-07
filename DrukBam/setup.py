#import pathlib
from setuptools import setup

# The directory containing this file
#HERE = pathlib.Path(__file__).parent

# The text of the README file
#README = (HERE / "README.md").read_text()

# This call to setup() does all the work
setup(
    name="drukbam",
    version="0.1",
    description="plots sorted indexed bam files",
    long_description='test',
    long_description_content_type="text/markdown",
    url="https://github.com/StephanHolgerD/DrukBam",
    author="Stephan Holger Drukewitz",
    author_email="steph-druk@web.de",
    license="MIT",
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
    ],
    packages=["DrukBamon"],
    include_package_data=True,
    install_requires=["pandas", "matplotlib",'pysam','multiprocessing','tqdm',],
    entry_points={
        "console_scripts": [
            "DrukBam=DrukBam.__main__:main",
        ]
    },
)
