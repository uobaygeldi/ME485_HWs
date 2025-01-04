# -*- coding: utf-8 -*-
from utils.nb import dot
from utils.np import eps

import numpy as np
import re


def get_fsolver(name, be, cplargs):
    """
    docstring
    """
    fname = re.sub(r'\+', 'p', name)
    fname = re.sub('-', '_', fname)
    flux = eval('make_' + fname)(cplargs)
    return be.compile(flux)


def make_rusanov(cplargs):
    nvars     = cplargs['nfvars']
    gamma     = cplargs['gamma']
    flux_func = cplargs['flux_func']
    array     = cplargs['array']
    ndims     = cplargs['ndims']
    def fsolver(ul, ur, vl, vr, nf, fn):
        fl = array(nvars)
        fr = array(nvars)

        # this is u*phi*n 
        flux_func(ul, vl, nf, fl)
        flux_func(ur, vr, nf, fr)
        for i in range(nvars):
            alpha = max(ul[i], ur[i])
            fn[i] = 0.5 * (fl[i] + fr[i]) - abs(alpha) * (ur[i] - ul[i])
        #---------------------------------#  
        # complete the function
        #---------------------------------#  

    return fsolver

def make_upwind(cplargs):
    nvars     = cplargs['nfvars']
    gamma     = cplargs['gamma']
    flux_func = cplargs['flux_func']
    array     = cplargs['array']
    ndims     = cplargs['ndims']
    def fsolver(ul, ur, vl, vr, nf, fn):
        fl = array(nvars)
        fr = array(nvars)

        # this is u*phi*n
        flux_func(ul, vl, nf, fl)
        flux_func(ur, vr, nf, fr)

        if dot(vl, nf, ndims) > 0:
            for i in range(nvars):
                fn[i] = fl[i]
        else:
            for i in range(nvars):
                fn[i] = fr[i]

    return fsolver