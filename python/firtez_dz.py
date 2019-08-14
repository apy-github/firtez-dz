#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

from scipy.io import FortranFile
import matplotlib.pyplot as pl

################################################################################
# Stokes profiles:
################################################################################

def read_profile(fname, fmt_type=np.float32):

  f = FortranFile(fname, 'r')
  first_rec = f.read_record(dtype=fmt_type)

  posv = first_rec[0]
  negv = first_rec[1]
  fid = first_rec[2]
  nrec = first_rec[3]
  ndims = first_rec[4]

  medv = (posv + negv) / 2.
  verp = posv - medv
  vern = negv - medv
  vers = (posv - negv) / 2.

  if ( (np.abs(medv-3000.) > 1.e-3) | (np.abs(fid-160904.) > 1.e-3) ):
    print('Something is wrong with the file %s' % fname)
    print('\tIs it a 3D profile file?')
    return np.nan

  if (np.abs(vers-3.) < 1.e-3):
    #VERSION 3
    profile3D = read_profile3D_v3(f,nrec,first_rec[5:],fmt_type)
  else:
    print('Version %i of profile file is not supported' % np.int(vers))
    return np.nan

  f.close()

  return profile3D

def read_profile3D_v3(ofile, nrec, dims, fmt_type):

  from pdb import set_trace as stop

  nx, ny, nw, ns = np.int16(dims)

  profile3D = stk_profile3D(nw, nx, ny)

  to_store = np.zeros(np.int64((ns*1.)*nx*ny*nw), dtype=np.float32)

  offset=0
  for it_nnn in range(np.int16(nrec)-1):
    tmp = ofile.read_record(dtype=fmt_type)
    #print tmp.size
    if (it_nnn == 0):
      indx = tmp * 1.
    if (it_nnn == 1):
      wave = tmp * 1.
    if (it_nnn >= 2):
      to_store[offset:offset+tmp.size] = tmp*1.
      offset=offset+tmp.size

  toload = np.moveaxis(to_store.reshape(nx,ny,nw,4)\
      , [0,1,2,3], [2,3,1,0]) * 1.
  #print toload.shape, (ns,nw,nx,ny), indx, wave
  profile3D.set_profiles(toload, indx, wave)

  return profile3D

class stk_profile3D(object):

  def __init__(self, nw, nx, ny):

    self.nw = nw * 1
    self.nx = nx * 1
    self.ny = ny * 1
    self.indx = np.zeros(self.nw)
    self.wave = np.zeros(self.nw)
    self.stki = np.zeros((self.nw, self.nx, self.ny))
    self.stkq = np.zeros((self.nw, self.nx, self.ny))
    self.stku = np.zeros((self.nw, self.nx, self.ny))
    self.stkv = np.zeros((self.nw, self.nx, self.ny))

    return

  def set_indx(self,array):
    self.indx = self.indx * 0. + array * 1.
    return

  def set_wave(self,array):
    self.wave = self.wave * 0. + array * 1.
    return

  def set_stki(self,array):
    self.stki = self.stki * 0. + array * 1.
    return

  def set_stkq(self,array):
    self.stkq = self.stkq * 0. + array * 1.
    return

  def set_stku(self,array):
    self.stku = self.stku * 0. + array * 1.
    return

  def set_stkv(self,array):
    self.stkv = self.stkv * 0. + array * 1.
    return

  def set_profiles(self, array, indx, wave):

    dims = np.array(array.shape)
    if ( (dims.size != 4) | (dims[1] != self.nw) | (dims[2] != self.nx)\
        | (dims[3] != self.ny) | (indx.size != self.nw) | (wave.size != self.nw) ):
      print('Wrong format')
      return

    self.set_indx(indx)
    self.set_wave(wave)
    self.set_stki(array[0,:,:,:])
    self.set_stkq(array[1,:,:,:])
    self.set_stku(array[2,:,:,:])
    self.set_stkv(array[3,:,:,:])

    return

  def set_to_write(self,v=None):

    if (v==None):
      array = np.zeros((4, self.nw, self.ny, self.nx))
      array[0,:,:,:] = np.moveaxis(self.stki, [0,1,2], [2,0,1]).T * 1.
      array[1,:,:,:] = np.moveaxis(self.stkq, [0,1,2], [2,0,1]).T * 1.
      array[2,:,:,:] = np.moveaxis(self.stku, [0,1,2], [2,0,1]).T * 1.
      array[3,:,:,:] = np.moveaxis(self.stkv, [0,1,2], [2,0,1]).T * 1.
    else:
      array = np.zeros((self.nx, self.ny, self.nw, 4))
      array[:,:,:,0] = np.moveaxis(self.stki, [0,1,2], [2,0,1]) * 1.
      array[:,:,:,1] = np.moveaxis(self.stkq, [0,1,2], [2,0,1]) * 1.
      array[:,:,:,2] = np.moveaxis(self.stku, [0,1,2], [2,0,1]) * 1.
      array[:,:,:,3] = np.moveaxis(self.stkv, [0,1,2], [2,0,1]) * 1.

    return self.indx, self.wave, array

  def write_profile(self, fname, vv=3, fmt_type=np.float32):

    if (vv == 3):
      write_profile3D_v3(fname, self.set_to_write(v=3), fmt_type=fmt_type)
    else:
      print('Version %i for profile file not supported!' % i (vv, ))

    return

  def plot(self,fnum=1,itx=[0,],ity=[0,], axis='w' \
      , pkwargs={}, fkwargs={}):

    if (type(itx)!=list):
      itx=list(itx)
    if (type(ity)!=list):
      ity=list(ity)

    if (axis=='p'):
      xtoplot = np.arange(self.nw)
    elif(axis=='w'):
      xtoplot = self.wave * 1.
    else:
      print("axis='p' or axis='w'")
      return

    fg,ax=pl.subplots(ncols=2,nrows=2,num=fnum,sharex=True, **fkwargs)
    for it_nnn in range(len(itx)):
      ax[0,0].plot(xtoplot, self.stki[:,itx[it_nnn],ity[it_nnn]], **pkwargs)
      ax[1,0].plot(xtoplot, self.stkq[:,itx[it_nnn],ity[it_nnn]], **pkwargs)
      ax[0,1].plot(xtoplot, self.stku[:,itx[it_nnn],ity[it_nnn]], **pkwargs)
      ax[1,1].plot(xtoplot, self.stkv[:,itx[it_nnn],ity[it_nnn]], **pkwargs)
    #ax[1,0].xaxis.set_ticks_position('right')

    if (np.mean(self.stki)<5.):
      ax[0,0].set_ylabel(r'I/I$_{{\rm c}}$')
      ax[0,1].set_ylabel(r'Q/I$_{{\rm c}}$')
      ax[1,0].set_ylabel(r'U/I$_{{\rm c}}$')
      ax[1,1].set_ylabel(r'V/I$_{{\rm c}}$')
    else:
      ax[0,0].set_ylabel(r'I')
      ax[0,1].set_ylabel(r'Q')
      ax[1,0].set_ylabel(r'U')
      ax[1,1].set_ylabel(r'V')

    if (axis == 'p'):
      ax[1,0].set_xlabel(r'$\lambda$ [px]')
      ax[1,1].set_xlabel(r'$\lambda$ [px]')

    elif (axis == 'w'):
      ax[1,0].set_xlabel(r'$\lambda$ [m$\AA$]')
      ax[1,1].set_xlabel(r'$\lambda$ [m$\AA$]')

      fg.tight_layout()

    return

def plot_profiles(profiles,fnum=1,itx=[0,],ity=[0,], axis='w', labels=[] \
    , linestyle=[], rangex=[] \
    , pargs=(), pkargs={}, fkwargs={}, pkwargs={}):

  if (type(itx)!=list):
    itx=list(itx)
  if (type(ity)!=list):
    ity=list(ity)
  if (type(profiles)!=list):
    profiles=list(profiles)
  if (type(labels)!=list):
    labels=list(labels)
  if (type(linestyle)!=list):
    linestyle=list(linestyle)

  if ( (len(rangex)!=0) & (len(rangex)!=2) ):
    print('rangex must be a list of two elements or not supplied')
    return

  show_labels = True
  if (len(labels)==0):
    show_labels = False
    labels = [''] * len(profiles)

  if (len(linestyle)==1):
    linestyle = linestyle * len(profiles)
  if (len(linestyle)==0):
    linestyle = ['-'] * len(profiles)

  if ( (len(labels)!=len(profiles)) ):
    print("Length of labels must be equal to profiles' length")
    print("\tLabels length: %i" % (len(labels),))
    print("\tLabels profiles: %i" % (len(profiles),))
    return

  fg,ax=pl.subplots(ncols=2,nrows=2,num=fnum,sharex=True, **fkwargs)
  lrangex=1.e99
  urangex=-1.e99
  for itn, self in enumerate(profiles):

    if (axis=='p'):
      xtoplot = np.arange(self.nw)
    elif(axis=='w'):
      xtoplot = self.wave * 1.
    else:
      print("axis='p' or axis='w'")
      return

    if (np.min(xtoplot)<lrangex):
      lrangex = np.min(xtoplot)
    if (np.max(xtoplot)>urangex):
      urangex = np.max(xtoplot)
  
    for it_nnn in range(len(itx)):

      ax[0,0].plot(xtoplot \
          , self.stki[:,itx[it_nnn],ity[it_nnn]], **pkwargs)
      ax[1,0].plot(xtoplot \
          , self.stkq[:,itx[it_nnn],ity[it_nnn]], **pkwargs)
      ax[0,1].plot(xtoplot \
          , self.stku[:,itx[it_nnn],ity[it_nnn]], **pkwargs)
      ax[1,1].plot(xtoplot \
          , self.stkv[:,itx[it_nnn],ity[it_nnn]], label=labels[itn] \
          , **pkwargs)

  if (np.mean(profiles[0].stki)<5.):
    ax[0,0].set_ylabel(r'I/I$_{{\rm c}}$')
    ax[0,1].set_ylabel(r'Q/I$_{{\rm c}}$')
    ax[1,0].set_ylabel(r'U/I$_{{\rm c}}$')
    ax[1,1].set_ylabel(r'V/I$_{{\rm c}}$')
  else:
    ax[0,0].set_ylabel(r'I')
    ax[0,1].set_ylabel(r'Q')
    ax[1,0].set_ylabel(r'U')
    ax[1,1].set_ylabel(r'V')

  if (axis == 'p'):
    ax[1,0].set_xlabel(r'$\lambda$ [px]')
    ax[1,1].set_xlabel(r'$\lambda$ [px]')

  elif (axis == 'w'):
    ax[1,0].set_xlabel(r'$\lambda$ [m$\AA$]')
    ax[1,1].set_xlabel(r'$\lambda$ [m$\AA$]')

  if (len(rangex)==2):
    lrangex=rangex[0]*1.
    urangex=rangex[1]*1.

    for it_xxx in range(3):
      for it_yyy in range(3):
        ax[it_xxx, it_yyy].set_xlim(lrangex, urangex)

  fg.tight_layout()
  if (show_labels):
    ax[1,1].legend(shadow=True, fancybox=True \
        , ncol=np.int(np.round(np.sqrt(len(profiles)))))
#, bbox_to_anchor=(0.5,-1.) \

  return





################################################################################
# Model:
################################################################################

def read_model(fname, fmt_type=np.float32):

  f = FortranFile(fname, 'r')

  first_rec = f.read_record(dtype=fmt_type)

  posv = first_rec[0]
  negv = first_rec[1]
  fid = first_rec[2]
  nrec = first_rec[3]
  ndims = first_rec[4]

  medv = (posv + negv) / 2.
  verp = posv - medv
  vern = negv - medv
  vers = (posv - negv) / 2.

  if ( (np.abs(medv-3000.) > 1.e-3) | (np.abs(fid-130904.) > 1.e-3) ):
    print('Something is wrong with the file %s' % fname)
    print('\tIs it a 3D model file?')
    return np.nan

  if (np.abs(vers-3.) < 1.e-3):
    #VERSION 3
    model3D = read_model3D_v3(f,nrec,first_rec[5:],fmt_type)
  else:
    print('Version %i for model atmosphere file not supported!' % np.int(vers))
    return np.nan

  f.close()

  return model3D

def read_model3D_v3(ofile, nrec, dims, fmt_type, verbose=False):

  npar, nx, ny, nz = np.int64(dims)

  if (np.abs(npar-13.) > 1.e-3):
    print('Something is wrong with the file %s' % fname)
    print('\tIs it a 3D model file?')
    return np.nan

  model3D = atm_model3D(nx, ny, nz)
  to_store = np.zeros(np.int((npar*1.)*nx*ny*nz), dtype=np.float32)

  offset=0
  for it_nnn in range(np.int64(nrec)-1):
    tmp = ofile.read_record(dtype=fmt_type)
    to_store[offset:offset+tmp.size] = tmp*1.
    offset=offset+tmp.size

  model3D.set_atmosphere(to_store.reshape(npar,nx,ny,nz))
  if (verbose==True):
    print('reading model: %s' % (fname, ))

  return model3D

class atm_model3D(object):

  def __init__(self, nx, ny, nz):

    self.nx = nx * 1
    self.ny = ny * 1
    self.nz = nz * 1
    self.tem = np.zeros((self.nx, self.ny, self.nz))
    self.pg = np.zeros((self.nx, self.ny, self.nz))
    self.rho = np.zeros((self.nx, self.ny, self.nz))
    self.bx = np.zeros((self.nx, self.ny, self.nz))
    self.by = np.zeros((self.nx, self.ny, self.nz))
    self.bz = np.zeros((self.nx, self.ny, self.nz))
    self.vz = np.zeros((self.nx, self.ny, self.nz))
    self.pel = np.zeros((self.nx, self.ny, self.nz))
    self.mw = np.zeros((self.nx, self.ny, self.nz))
    self.tau = np.zeros((self.nx, self.ny, self.nz))
    self.x = np.zeros((self.nx, self.ny, self.nz))
    self.y = np.zeros((self.nx, self.ny, self.nz))
    self.z = np.zeros((self.nx, self.ny, self.nz))

    return

  def set_tem(self,array):
    self.tem = self.tem * 0. + array * 1.
    return

  def set_pg(self,array):
    self.pg = self.pg * 0. + array * 1.
    return

  def set_rho(self,array):
    self.rho = self.rho * 0. + array * 1.
    return

  def set_bx(self,array):
    self.bx = self.bx * 0. + array * 1.
    return

  def set_by(self,array):
    self.by = self.by * 0. + array * 1.
    return

  def set_bz(self,array):
    self.bz = self.bz * 0. + array * 1.
    return

  def set_vz(self,array):
    self.vz = self.vz * 0. + array * 1.
    return

  def set_pel(self,array):
    self.pel = self.pel * 0. + array * 1.
    return

  def set_mw(self,array):
    self.mw = self.mw * 0. + array * 1.
    return

  def set_tau(self,array):
    self.tau = self.tau * 0. + array * 1.
    return

  def set_x(self,array):
    self.x = self.x * 0. + array * 1.
    return

  def set_y(self,array):
    self.y = self.y * 0. + array * 1.
    return

  def set_z(self,array):
    self.z = self.z * 0. + array * 1.
    return

  def set_atmosphere(self, array):

    dims = np.array(array.shape)

    if ( (dims.size != 4)\
        | (dims[1] != self.nx)\
        | (dims[2] != self.ny)\
        | (dims[3] != self.nz)\
        ):
      print('Wrong format')
      return

    self.set_tem(array[0,:,:,:])
    self.set_pg(array[1,:,:,:])
    self.set_rho(array[2,:,:,:])
    self.set_bx(array[3,:,:,:])
    self.set_by(array[4,:,:,:])
    self.set_bz(array[5,:,:,:])
    self.set_vz(array[6,:,:,:])
    self.set_pel(array[7,:,:,:])
    self.set_mw(array[8,:,:,:])
    self.set_tau(array[9,:,:,:])
    self.set_x(array[10,:,:,:])
    self.set_y(array[11,:,:,:])
    self.set_z(array[12,:,:,:])

    return

  def set_to_write(self,v=None):

    if (v==None):
      array = np.zeros((13, self.nz, self.ny, self.nx))

      array[0,:,:,:] = self.tem.T * 1.
      array[1,:,:,:] = self.pg.T * 1.
      array[2,:,:,:] = self.rho.T * 1.
      array[3,:,:,:] = self.bx.T * 1.
      array[4,:,:,:] = self.by.T * 1.
      array[5,:,:,:] = self.bz.T * 1.
      array[6,:,:,:] = self.vz.T * 1.
      array[7,:,:,:] = self.pel.T * 1.
      array[8,:,:,:] = self.mw.T * 1.
      array[9,:,:,:] = self.tau.T * 1.
      array[10,:,:,:] = self.x.T * 1.
      array[11,:,:,:] = self.y.T * 1.
      array[12,:,:,:] = self.z.T * 1.
    elif(v==3):
      array = np.zeros((13, self.nx, self.ny, self.nz))

      array[0,:,:,:] = self.tem * 1.
      array[1,:,:,:] = self.pg * 1.
      array[2,:,:,:] = self.rho * 1.
      array[3,:,:,:] = self.bx * 1.
      array[4,:,:,:] = self.by * 1.
      array[5,:,:,:] = self.bz * 1.
      array[6,:,:,:] = self.vz * 1.
      array[7,:,:,:] = self.pel * 1.
      array[8,:,:,:] = self.mw * 1.
      array[9,:,:,:] = self.tau * 1.
      array[10,:,:,:] = self.x * 1.
      array[11,:,:,:] = self.y * 1.
      array[12,:,:,:] = self.z * 1.

    return array

  def write_model(self, fname, vv=3, fmt_type=np.float32, verbose=False):

    if (verbose==True):
      print('Writing model: %s' % (fname,))
    if (vv == 3):
      write_model3D_v3(fname, self.set_to_write(v=vv), fmt_type=fmt_type)
    else:
      print('Version %i for model atmosphere not supported!' % (vv,))
      return np.nan
    return


  def re_sample(self, new_nz, step=0.\
      , kind='linear', fill_value='extrapolate', new_zz=np.array([])):

    if (np.abs(step) > 0.1):
      #.step = (self.z[-1] - self.z[0]) / np.float32(new_nz)
      new_nz = np.int16(np.floor((self.z[-1] - self.z[0]) / np.float32(step)))

    old_z = self.z[0, 0, :] * 1.
    new_z = np.linspace(self.z[0,0,0], self.z[0,0,-1], new_nz)

    if (new_zz.size != 0):
      new_nz = new_zz.size
      new_z = new_zz * 1.

    new_model = atm_model3D(self.nx, self.ny, new_nz)

    narray = np.zeros((13, self.nx, self.ny, new_nz))
    for it_nx in range(self.nx):
      for it_ny in range(self.ny):

        f = interp1d(old_z, self.tem[it_nx, it_ny, :]\
            , kind=kind, fill_value=fill_value)
        narray[0, it_nx, it_ny, :] = f(new_z)
        f = interp1d(old_z, self.pg[it_nx, it_ny, :]\
            , kind=kind, fill_value=fill_value)
        narray[1, it_nx, it_ny, :] = f(new_z)
        f = interp1d(old_z, self.rho[it_nx, it_ny, :]\
            , kind=kind, fill_value=fill_value)
        narray[2, it_nx, it_ny, :] = f(new_z)
        f = interp1d(old_z, self.bx[it_nx, it_ny, :]\
            , kind=kind, fill_value=fill_value)
        narray[3, it_nx, it_ny, :] = f(new_z)
        f = interp1d(old_z, self.by[it_nx, it_ny, :]\
            , kind=kind, fill_value=fill_value)
        narray[4, it_nx, it_ny, :] = f(new_z)
        f = interp1d(old_z, self.bz[it_nx, it_ny, :]\
            , kind=kind, fill_value=fill_value)
        narray[5, it_nx, it_ny, :] = f(new_z)
        f = interp1d(old_z, self.vz[it_nx, it_ny, :]\
            , kind=kind, fill_value=fill_value)
        narray[6, it_nx, it_ny, :] = f(new_z)
        f = interp1d(old_z, self.pel[it_nx, it_ny, :]\
            , kind=kind, fill_value=fill_value)
        narray[7, it_nx, it_ny, :] = f(new_z)
        f = interp1d(old_z, self.mw[it_nx, it_ny, :]\
            , kind=kind, fill_value=fill_value)
        narray[8, it_nx, it_ny, :] = f(new_z)
        f = interp1d(old_z, self.tau[it_nx, it_ny, :]\
            , kind=kind, fill_value=fill_value)
        narray[9, it_nx, it_ny, :] = f(new_z)
        f = interp1d(old_z, self.x[it_nx, it_ny, :]\
            , kind=kind, fill_value=fill_value)
        narray[10, it_nx, it_ny, :] = f(new_z)
        f = interp1d(old_z, self.y[it_nx, it_ny, :]\
            , kind=kind, fill_value=fill_value)
        narray[11, it_nx, it_ny, :] = f(new_z)
        narray[12, it_nx, it_ny, :] = new_z * 1.

    new_model.set_atmosphere(narray)

    return new_model

  def plot(self,fnum=1,itx=[0,],ity=[0,], axis='z', pargs=(), pkargs={}, fkwargs={}):

    if (type(itx)!=list):
      itx=list(itx)
    if (type(ity)!=list):
      ity=list(ity)

    fg,ax=pl.subplots(ncols=3,nrows=3,num=fnum, **fkwargs)

    if (axis=='z'):
      xtoplot = self.z * 1.e-3
    elif(axis=='t'):
      xtoplot = self.tau * 1.
    else:
      print("axis='z' or axis='t'")
      return

    for it_nnn in range(len(itx)):
      ax[0,0].step(xtoplot[itx[it_nnn],ity[it_nnn],:], self.tem[itx[it_nnn],ity[it_nnn],:] * 1.e-3)
      ax[0,1].step(xtoplot[itx[it_nnn],ity[it_nnn],:], self.pg[itx[it_nnn],ity[it_nnn],:])
      ax[0,1].set_yscale('log')
      ax[0,2].step(xtoplot[itx[it_nnn],ity[it_nnn],:], self.rho[itx[it_nnn],ity[it_nnn],:])
      ax[0,2].set_yscale('log')
      ax[1,0].step(xtoplot[itx[it_nnn],ity[it_nnn],:], self.bx[itx[it_nnn],ity[it_nnn],:] * 1.e-3)
      ax[1,1].step(xtoplot[itx[it_nnn],ity[it_nnn],:], self.by[itx[it_nnn],ity[it_nnn],:] * 1.e-3)
      ax[1,2].step(xtoplot[itx[it_nnn],ity[it_nnn],:], self.bz[itx[it_nnn],ity[it_nnn],:] * 1.e-3)
      ax[2,0].step(xtoplot[itx[it_nnn],ity[it_nnn],:], self.vz[itx[it_nnn],ity[it_nnn],:]*1.e-5)
      if (axis == 'z'):
        ax[2,1].step(xtoplot[itx[it_nnn],ity[it_nnn],:], self.tau[itx[it_nnn],ity[it_nnn],:])
      elif (axis=='t'):
        ax[2,1].step(xtoplot[itx[it_nnn],ity[it_nnn],:], self.z[itx[it_nnn],ity[it_nnn],:])
    fg.delaxes(ax[2,2])


    ax[0,0].xaxis.set_ticklabels('')
    ax[0,1].xaxis.set_ticklabels('')
    ax[0,2].xaxis.set_ticklabels('')
    ax[1,0].xaxis.set_ticklabels('')
    ax[1,1].xaxis.set_ticklabels('')

    ax[0,0].set_ylabel(r'T [Kk]')
    ax[0,1].set_ylabel(r'P$_{g}$ [dyn/cm$^2$]')
    ax[0,2].set_ylabel(r'$\rho$ [gr/$cm^3$]')
    ax[1,0].set_ylabel(r'B$_{x}$ [KG]')
    ax[1,1].set_ylabel(r'B$_{y}$ [KG]')
    ax[1,2].set_ylabel(r'B$_{z}$ [KG]')
    ax[2,0].set_ylabel(r'v$_{z}$ [km/s]')
    if (axis == 'z'):
      ax[2,1].set_ylabel(r'$\lg\tau_{5}$')

      ax[1,2].set_xlabel(r'z [Mm]')
      ax[2,0].set_xlabel(r'z [Mm]')
      ax[2,1].set_xlabel(r'z [Mm]')
    elif (axis == 't'):
      ax[2,1].set_ylabel(r'z [Mm]')

      ax[1,2].set_xlabel(r'$\lg\tau_{5}$')
      ax[2,0].set_xlabel(r'$\lg\tau_{5}$')
      ax[2,1].set_xlabel(r'$\lg\tau_{5}$')

    fg.tight_layout()

    return

def write_model3D_v3(fname, data, fmt_type=np.float32):

  npar, nx, ny, nz = data.shape
  #print('Writing model:')
  #print(data.reshape(npar,-1).sum(1)/np.float(nx)/np.float(ny)/np.float(nz))
  data = data.flatten()

  idd = 130904
  v=3

  if (fmt_type != np.float32):
    print('Not implemented yet!')
    return np.nan

  # We write single precision numbers: 32 bits
  # This is 4 bytes

  nnum_to_write = (npar * 1.) * nz * ny * nx
  nbyts_to_write = nnum_to_write * 4.

  #print(nnum_to_write)

  # The maximum size of a fortran record is 2Gb
  nbyts_max_frec = 2. * 1024. * 1024. * 1024.
  nnum_max_frec = nbyts_max_frec / 4.

  # We set the maximum number of fortran rec to actual_max-10:
  nnum_rec = np.int(nnum_max_frec - 10.)
  nrecs_to_write = np.int(np.ceil(nnum_to_write/nnum_rec))

  # FIRST, WE WRITE A FIRST RECORD WITH THE MANDATORY DATA:
  towrite = np.zeros(9, dtype=fmt_type)

  towrite[0] = 3000 + v
  towrite[1] = 3000 - v
  towrite[2] = idd * 1.
  towrite[3] = nrecs_to_write + 1. #NREC
  towrite[4] = 4.
  towrite[5] = npar * 1.
  towrite[6] = nx * 1.
  towrite[7] = ny * 1.
  towrite[8] = nz * 1.

  #OPEN WRITING FILE:
  f=FortranFile(fname,'w')
  # FIRST ACCESS:
  f.write_record(np.float32(towrite))
  # ONGOING ACCESES:
  #         SET OF DATA:
  offset = 0
  for it_nnn in range(nrecs_to_write):
    #print('Writing %i of %i' % (it_nnn+1, nrecs_to_write))
    if (it_nnn != nrecs_to_write - 1):
      towrite = data[offset:offset+nnum_rec] * 1.
      #print towrite.size, nnum_rec
    else:
      towrite = data[offset:] * 1.
      #print towrite.size, nnum_rec
    f.write_record(np.float32(towrite))
    offset=offset+nnum_rec
  # CLOSE IT
  f.close()

  return
#
def plot_models(models,fnum=1,itx=[0,],ity=[0,], axis='z', labels=[] \
    , linestyle=[], rangex=[] \
    , pargs=(), pkargs={}, fkwargs={}):

  if (type(itx)!=list):
    itx=list(itx)
  if (type(ity)!=list):
    ity=list(ity)
  if (type(models)!=list):
    models=list(models)
  if (type(labels)!=list):
    labels=list(labels)
  if (type(linestyle)!=list):
    linestyle=list(linestyle)

  if ( (len(rangex)!=0) & (len(rangex)!=2) ):
    print('rangex must be a list of two elements or not supplied')
    return

  show_labels = True
  if (len(labels)==0):
    show_labels = False
    labels = [''] * len(models)

  if (len(linestyle)==1):
    linestyle = linestyle * len(models)
  if (len(linestyle)==0):
    linestyle = ['-'] * len(models)

  if ( (len(labels)!=len(models)) ):
    print("Length of labels must be equal to models' length")
    print("\tLabels length: %i" % (len(labels),))
    print("\tLabels models: %i" % (len(models),))
    return

  fg,ax=pl.subplots(ncols=3,nrows=3,num=fnum, **fkwargs)
  lrangex=1.e99
  urangex=-1.e99
  for itn, self in enumerate(models):

    if (axis=='z'):
      xtoplot = self.z * 1.e-3
    elif(axis=='t'):
      xtoplot = self.tau * 1.
    else:
      print("axis='z' or axis='t'")
      return

    if (np.min(xtoplot)<lrangex):
      lrangex = np.min(xtoplot)
    if (np.max(xtoplot)>urangex):
      urangex = np.max(xtoplot)
  
    for it_nnn in range(len(itx)):
      ax[0,0].step(xtoplot[itx[it_nnn],ity[it_nnn],:] \
          , self.tem[itx[it_nnn],ity[it_nnn],:] * 1.e-3 \
          , linestyle=linestyle[itn])
      ax[0,1].step(xtoplot[itx[it_nnn],ity[it_nnn],:] \
          , self.pg[itx[it_nnn],ity[it_nnn],:] \
          , linestyle=linestyle[itn])
      ax[0,1].set_yscale('log')
      ax[0,2].step(xtoplot[itx[it_nnn],ity[it_nnn],:] \
          , self.rho[itx[it_nnn],ity[it_nnn],:] \
          , linestyle=linestyle[itn])
      ax[0,2].set_yscale('log')
      ax[1,0].step(xtoplot[itx[it_nnn],ity[it_nnn],:] \
          , self.bx[itx[it_nnn],ity[it_nnn],:] * 1.e-3 \
          , linestyle=linestyle[itn])
      ax[1,1].step(xtoplot[itx[it_nnn],ity[it_nnn],:] \
          , self.by[itx[it_nnn],ity[it_nnn],:] * 1.e-3 \
          , linestyle=linestyle[itn])
      ax[1,2].step(xtoplot[itx[it_nnn],ity[it_nnn],:] \
          , self.bz[itx[it_nnn],ity[it_nnn],:] * 1.e-3, label=labels[itn] \
          , linestyle=linestyle[itn])
      ax[2,0].step(xtoplot[itx[it_nnn],ity[it_nnn],:] \
          , self.vz[itx[it_nnn],ity[it_nnn],:]*1.e-5 \
          , linestyle=linestyle[itn])
      if (axis == 'z'):
        ax[2,1].step(xtoplot[itx[it_nnn],ity[it_nnn],:] \
          , self.tau[itx[it_nnn],ity[it_nnn],:] \
          , linestyle=linestyle[itn])
      elif (axis=='t'):
        ax[2,1].step(xtoplot[itx[it_nnn],ity[it_nnn],:] \
          , self.z[itx[it_nnn],ity[it_nnn],:] \
          , linestyle=linestyle[itn])

  ax[0,1].xaxis.set_ticklabels('')
  ax[0,2].xaxis.set_ticklabels('')
  ax[1,0].xaxis.set_ticklabels('')
  ax[1,1].xaxis.set_ticklabels('')

  ax[0,0].set_ylabel(r'T [Kk]')
  ax[0,1].set_ylabel(r'P$_{g}$ [dyn/cm$^2$]')
  ax[0,2].set_ylabel(r'$\rho$ [gr/$cm^3$]')
  ax[1,0].set_ylabel(r'B$_{x}$ [KG]')
  ax[1,1].set_ylabel(r'B$_{y}$ [KG]')
  ax[1,2].set_ylabel(r'B$_{z}$ [KG]')
  ax[2,0].set_ylabel(r'v$_{z}$ [km/s]')
  if (axis == 'z'):
    ax[2,1].set_ylabel(r'$\lg\tau_{5}$')

    ax[1,2].set_xlabel(r'z [Mm]')
    ax[2,0].set_xlabel(r'z [Mm]')
    ax[2,1].set_xlabel(r'z [Mm]')
  elif (axis == 't'):
    ax[2,1].set_ylabel(r'z [Mm]')

    ax[1,2].set_xlabel(r'$\lg\tau_{5}$')
    ax[2,0].set_xlabel(r'$\lg\tau_{5}$')
    ax[2,1].set_xlabel(r'$\lg\tau_{5}$')

  if (len(rangex)==2):
    lrangex=rangex[0]*1.
    urangex=rangex[1]*1.

    for it_xxx in range(3):
      for it_yyy in range(3):
        ax[it_xxx, it_yyy].set_xlim(lrangex, urangex)
  fg.delaxes(ax[2,2])

  ax[0,0].xaxis.set_ticklabels('')

  fg.tight_layout()
  if (show_labels):
    ax[1,2].legend(loc=8,shadow=True, fancybox=True, bbox_to_anchor=(0.5,-1.) \
        , ncol=np.int(np.round(np.sqrt(len(models)))))

  return



