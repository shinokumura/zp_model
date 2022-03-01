#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import codecs
from setuptools import setup


def read(fname):
    file_path = os.path.join(os.path.dirname(__file__), fname)
    return codecs.open(file_path, encoding="utf-8").read()


setup(
    name="pytest-zpmodel",
    version="0.1.0",
    author="Shin Okumura",
    author_email="s.okumura@iaea.org",
    maintainer="Shin Okumura",
    maintainer_email="s.okumura@iaea.org",
    license="MIT",
    url="https://github.com/s.okumura@iaea.org/pytest-zpmodel",
    description="Revised Zp model",
    long_description=read("README.rst"),
    py_modules=["pytest_zpmodel"],
    python_requires=">=3.5",
    install_requires=["pytest>=3.8.5"],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Framework :: Pytest",
        "Intended Audience :: Developers",
        "Topic :: Software Development :: Testing",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3 :: Only",
        "Programming Language :: Python :: Implementation :: CPython",
        "Programming Language :: Python :: Implementation :: PyPy",
        "Operating System :: OS Independent",
        "License :: OSI Approved :: MIT License",
    ],
    entry_points={"pytest11": ["zpmodel = pytest_zpmodel",],},
)
