#!/usr/bin/env python

# -*- coding: utf-8 -*-

"""mtl.__main__: executed when mtl directory is called as script."""

import sys
from mtl.mtl import run_by_environment_vars as main

if __name__ == '__main__':
    sys.exit(main())
