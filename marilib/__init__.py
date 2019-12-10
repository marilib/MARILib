#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 25 16:22:36 2018
.. moduleauthor:: François Gallard
.. moduleauthor:: Thales Delmiro
"""
import numpy as __np
from scipy.optimize import fsolve as __fs
import marilib
import sys

marilib.is_driven_by_gems = False

marilib.is_using_autograd = False
sys.modules["marilib"].numpy = __np
sys.modules["marilib.numpy"] = __np

marilib.is_using_newton_solve = False
sys.modules["marilib"].fsolve = __fs
sys.modules["marilib.fsolve"] = __fs

def use_newton_solve():
    from marilib.tools.math import newton_solve as __agfs
    sys.modules["marilib"].fsolve = __agfs
    sys.modules["marilib.fsolve"] = __agfs
    marilib.is_using_newton_solve = True

def use_autograd():
    from autograd import numpy as __agnp
    sys.modules["marilib"].numpy = __agnp
    sys.modules["marilib.numpy"] = __agnp
    marilib.is_using_autograd = True
    use_newton_solve()
