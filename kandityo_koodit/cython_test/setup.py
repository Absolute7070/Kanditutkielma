#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 08:31:54 2023

@author: degnaiyu
"""
from distutils.core import setup
from Cython.Build import cythonize
 
setup(
    ext_modules = cythonize("function.pyx")
)
