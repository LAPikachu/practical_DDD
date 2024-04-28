#!/usr/bin/env python

from numpy.distutils.core import setup,Extension

fsrc = Extension(name = 'vtkam.fsrc',
                 sources = ['vtkam/varglob.f90','vtkam/fortran_src.f90'],
                # f2py_options=['--debug-capi']
                )

setup (
name         = "VTKAM",
version      = "v2.0",
description  = '',
author       = 'LEM',
author_email = 'laurent.korzeczek@onera.fr',
packages     = ['vtkam'],
ext_modules = [fsrc]
)
