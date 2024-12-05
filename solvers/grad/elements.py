# -*- coding: utf-8 -*-
import math
from time import process_time
import numpy as np
import re

from solvers.base import BaseElements
from backends.types import ArrayBank, Kernel, NullKernel
from utils.np import eps, chop, npeval


class gradFluidElements:
    @property
    def primevars(self):
        # Primitive variables
        return ['q'] 
    @property
    def conservars(self):
        # Conservative variables
        pri = self.primevars

        # rho,rhou,rhov,rhow,E
        return [pri[0]]

    def prim_to_conv(self, pri, cfg):
        return [pri[0]]

    def conv_to_prim(self, con, cfg):
        return [con[0]]

    def fix_nonPys_container(self):
        # Constants and dimensions
        ndims, nfvars = self.ndims, self.nfvars

        def fix_nonPhy(u):
            u[0] = u[0]
        # Compile the function
        return self.be.compile(fix_nonPhy)

class GradElements(BaseElements,  gradFluidElements):
    nreg = 1
    def __init__(self, be, cfg, name, eles):
        super().__init__(be, cfg, name, eles)
        self.nvars = len(self.primevars)
        self.nfvars = self.nvars
        self._const = cfg.items('constants')
        self._grad_method = cfg.get('solver', 'gradient')
        # print(self.dxv)

    def construct_kernels(self, vertex, nreg, impl_op=0):
        self.vertex = vertex

        # Upts : Solution vector
        self.upts = upts = [self._ics.copy() for i in range(nreg)]
        del(self._ics)

        # Solution vector bank and assign upts index
        self.upts_in = upts_in = ArrayBank(upts, 0)
        self.upts_out = upts_out = ArrayBank(upts, 1)

        # Construct arrays for flux points, dt and derivatives of source term
        self.fpts = fpts = np.empty((self.nface, self.nvars, self.neles))
        self.grad = grad = np.zeros((self.ndims, self.nvars, self.neles))

        lim = np.ones((self.nvars, self.neles))

        limiter = self.cfg.get('solver', 'limiter', 'none')
        # Prepare vertex array
        vpts = vertex.make_array(limiter)
        
        self.compute_fpts = Kernel(self._make_compute_fpts(), upts_in, fpts)

       
        if( (self._grad_method=='least-square') or (self._grad_method=='weighted-least-square')):
            self.compute_grad = Kernel(self._make_grad_ls(), fpts, grad)       
        elif(self._grad_method == 'green-gauss-cell'):
            self.compute_grad = Kernel(self._make_grad_gg(), fpts, grad)
        else:
            self.compute_grad = Kernel(self._make_grad_ls(), fpts, grad)

        
        # Kernel for linear reconstruction
        self.compute_recon = Kernel(self._make_recon(), upts_in, grad, lim, fpts)
        
        # Kernel to post-process
        self.post = Kernel(self._make_post(), upts_in)
 #-------------------------------------------------------------------------------#
    def compute_L2_norm(self):
        # adjust faces to eles
        nvars, neles, ndims = self.nvars, self.neles, self.ndims
        vol                 = self._vol
        xcs                 = self.xc # Cell center point (x,y)
        # print(np.shape(xcs)) # (elem id, dimension)
        # print(np.shape(vol))
        # self.grad = grad = np.zeros((self.ndims, self.nvars, self.neles))
        exactGrad = np.zeros((self.ndims, self.nvars, self.neles))
        err = np.zeros((self.ndims, self.nvars, self.neles))
        # sums = np.zeros((self.ndims, self.nvars))
        for idx in range(neles):
            x = np.zeros(ndims)
            for i in range(nvars):
                for j in range(ndims):
                    x[j] = xcs[idx][j]
                eqn = [x[0]/math.sqrt(x[0]*x[0]+x[1]*x[1]), x[1]/math.sqrt(
                    x[0]*x[0] + x[1]*x[1])] # sqrt(x^2 + y^2)
                # eqn = [2*x[0],2*x[1]] # x^2 + y^2
                for j in range(ndims):
                    exactGrad[j, i, idx] = eqn[j]
                    err[j, i, idx] = (self.grad[j, i, idx] - exactGrad[j, i, idx]
                                       )*vol[idx]

        resid = np.linalg.norm(err, axis=2)
        return resid
#-------------------------------------------------------------------------------#
    # Assign cell centers values to face centers
    def _make_compute_fpts(self):
        nvars, nface = self.nvars, self.nface

        def _compute_fpts(i_begin, i_end, upts, fpts):
            # Copy upts to fpts
            # fpts = np.empty((self.nface, self.nvars, self.neles))
            # print(np.shape(upts))
            for idx in range(i_begin, i_end):
                for j in range(nvars):
                    for i in range(nface):
                        fpts[i, j, idx] = upts[j, idx]

        return self.be.make_loop(self.neles, _compute_fpts)   
#-------------------------------------------------------------------------------#
    def _make_grad_ls(self):
        nface, ndims, nvars = self.nface, self.ndims, self.nvars
        # Gradient operator 
        op = self._grad_operator
        # print(np.shape(op))
        def _cal_grad(i_begin, i_end, fpts, grad):
            # Elementwiseloop
            for idx in range(i_begin, i_end):
                for d in range(ndims):
                    for j in range(nvars):  # = 1
                        sum = 0
                        for i in range(nface):
                            # fpts = np.empty((self.nface, self.nvars, self.neles))
                            sum += fpts[i, j, idx] * op[d, i, idx]
                            #grad = np.zeros((self.ndims, self.nvars, self.neles))

                        grad[d, j, idx] = sum

        # Compile the function
        return self.be.make_loop(self.neles, _cal_grad)    
#-------------------------------------------------------------------------------#
    def _make_grad_gg(self):

        nface, ndims, nvars = self.nface, self.ndims, self.nvars
        # Normal vector and volume
        snorm_mag = self._mag_snorm_fpts # face | elemid
        snorm_vec = np.rollaxis(self._vec_snorm_fpts, 2) # x,y | face | elemid
        vol       = self._vol # elemid
        # fpts = np.empty((self.nface, self.nvars, self.neles))
        def _cal_grad(i_begin, i_end, fpts, grad):
            # Elementwise loop starts here
            # print(snorm_mag)
            # print(snorm_vec)
            for idx in range(i_begin, i_end):
                for d in range(ndims):
                    for j in range(nvars): # = 1
                        sum = 0
                        for i in range(nface):
                            sum += fpts[i, j, idx]*snorm_mag[i,idx]*snorm_vec[d,i,idx]
                            #grad = np.zeros((self.ndims, self.nvars, self.neles))

                        grad[d,j,idx] = (1/vol[idx])*sum

        return self.be.make_loop(self.neles, _cal_grad)      

#-------------------------------------------------------------------------------#
    def _make_recon(self):
        nface, ndims, nvars = self.nface, self.ndims, self.nvars

        # Displacement vector
        op = self.dxf

        def _cal_recon(i_begin, i_end, upts, grad, lim, fpts):
            # Elementwise dot product and scale with limiter
            # TODO: Reduce accesing global array
            for i in range(i_begin, i_end):
                for l in range(nvars):
                    for k in range(nface):
                        tmp = 0
                        for j in range(ndims):
                            tmp += op[k, j, i]*grad[j, l, i]
                                                    
                        fpts[k, l, i] = upts[l, i] + lim[l, i]*tmp

        return self.be.make_loop(self.neles, _cal_recon)

#-------------------------------------------------------------------------------#

    @property
    # @fc.lru_cache()
    # @chop
    
    def _grad_operator(self):
        # Difference of displacement vector (cell to cell)
        # (Nfaces, Nelements, dim) -> (dim, Nfaces, Nelements)
        dxc = np.rollaxis(self.dxc, 2)
        # (Nfaces, Nelements)
        distance = np.linalg.norm(dxc, axis=0)

        # Code taken directly from base/elements.py, thanks hocam :)

        # Normal vector and volume
        snorm_mag = self._mag_snorm_fpts
        snorm_vec = np.rollaxis(self._vec_snorm_fpts, 2)
        vol = self._vol

        if self._grad_method == 'least-square':
            beta, w = 1.0, 1.0
        elif self._grad_method == 'weighted-least-square':
            # Invserse distance weight
            beta, w = 1.0, 1 / distance ** 2
        elif self._grad_method == 'green-gauss':
            beta, w = 0.0, 1.0
        elif self._grad_method == 'hybrid':
            # Shima et al., Green–Gauss/Weighted-Least-Squares
            # Hybrid Gradient Reconstruction for
            # Arbitrary Polyhedra Unstructured Grids, AIAA J., 2013
            # WLSQ(G)
            dxf = self.dxf.swapaxes(0, 1)

            dxcn = np.einsum('ijk,ijk->jk', dxc, snorm_vec)
            dxfn = np.einsum('ijk,ijk->jk', dxf, snorm_vec)

            w = (2 * dxfn / dxcn) ** 2 * snorm_mag / distance

            # Compute blending function (GLSQ)
            ar = 2 * np.linalg.norm(self.dxf, axis=1).max(
                axis=0)*snorm_mag.max(axis=0) / vol
            beta = np.minimum(1, 2 / ar)
        else:
            raise ValueError("Invalid gradient method : ", self._grad_method)

        # Scaled dxc vector
        dxcs = dxc * np.sqrt(w)

        # Least square matrix [dx*dy] and its inverse
        lsq = np.array([[np.einsum('ij,ij->j', x, y)
                         for y in dxcs] for x in dxcs])

        # Hybrid type of Ax=b
        A = beta * lsq + 2 * (1 - beta) * vol * np.eye(self.ndims)[:, :, None]
        b = beta * dxc * w + 2 * (1 - beta) * 0.5 * snorm_vec * snorm_mag

        # Solve Ax=b
        op = np.linalg.solve(np.rollaxis(
            A, axis=2), np.rollaxis(b, axis=2)).transpose(1, 2, 0)

        return op

    def _make_post(self):
        nface, ndims, nvars = self.nface, self.ndims, self.nvars
        grad = self.grad
        xc = self.xc.T        
        def post(i_begin, i_end, upts):
            # Apply the function over eleemnts
            for idx in range(i_begin, i_end):          
                print(idx, grad[0, 0, idx], grad[1, 0, idx])

        return self.be.make_loop(self.neles, post)
