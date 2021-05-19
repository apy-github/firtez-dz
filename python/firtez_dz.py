#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#
# To remove:
#
from pdb import set_trace as stop
#
#
import numpy as np
#
from scipy.io import FortranFile
from scipy.interpolate import interp1d
import matplotlib.pyplot as pl
from matplotlib import ticker
pl.ion()
#
pl.rcParams.update({'font.size':15})
#pl.rcParams.update({'text.usetex':False})
pl.rcParams.update({'text.usetex':True})
pl.rcParams['xtick.minor.visible'] = True
pl.rcParams['ytick.minor.visible'] = True
#
import warnings
warnings.filterwarnings('ignore')
#
################################################################################
# Stokes profiles:
################################################################################

def read_profile(fname, fmt_type=np.float32):

  def read_profile3D_v3(ofile, nrec, dims, fmt_type):
  
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

class stk_profile3D(object):

  def __init__(self, nw, nx, ny):
    """
      Special method to initialize firtez-dz class:

        Inputs:
          1/3- nw: integer, number of wavelengths
          2/3- nx: integer, number of pixels in x direction
          3/3- ny: integer, number of pixels in y direction

    """

    self.nw = nw * 1
    self.nx = nx * 1
    self.ny = ny * 1
    self.indx = np.zeros(self.nw)
    self.wave = np.zeros(self.nw)
    self.stki = np.zeros((self.nw, self.nx, self.ny))
    self.stkq = np.zeros((self.nw, self.nx, self.ny))
    self.stku = np.zeros((self.nw, self.nx, self.ny))
    self.stkv = np.zeros((self.nw, self.nx, self.ny))
    self.shape = np.shape(self.stki)
    # For tv:
    self.__tv_defaults = {}
    self.__set_tv_defaults()

    return
  #
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  #
  def __repr__(self):
    """
      Special method 'repr':
    """

    return "\n\tFirtez-dz Stokes class."
  #
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  #
  def __add__(self, stkcls_to_add):
    """
      Special method :
    """

    #
    # Check alignment:
    if ( (self.nx!=stkcls_to_add.nx) | (self.ny!=stkcls_to_add.ny) ):
      print("")
      print("\tError!")
      print("\tProfiles to add are not aligned!")
      print("")
      return np.nan

    out = stk_profile3D(self.nw+stkcls_to_add.nw, self.nx, self.ny)

    for itp in ["indx","wave","stki","stkq","stku","stkv"]:
      if (not itp in out):
        #
        # Second 
        from pdb import set_trace as stop
        stop()
      out[itp] = np.concatenate([self[itp],stkcls_to_add[itp]], axis=0)
        

    return out
  #
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  #
  def __iadd__(self, stkcls_to_add):
    """
      Special method :
    """

    #
    # Check alignment:
    if ( (self.nx!=stkcls_to_add.nx) | (self.ny!=stkcls_to_add.ny) ):
      print("")
      print("\tError!")
      print("\tProfiles to add are not aligned!")
      print("")
      return np.nan

    for itp in ["indx","wave","stki","stkq","stku","stkv"]:
      if (not itp in self):
        #
        # Second 
        from pdb import set_trace as stop
        stop()
      self[itp] = np.concatenate([self[itp],stkcls_to_add[itp]], axis=0)
 
    self.shape = self.stki.shape
    self.nw = self.shape[0] * 1

    return self
  #
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  #
  def __getitem__(self, key):
    """
      Special method :
    """
    return getattr(self, key)
  #
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  #
  def __setitem__(self, key, item):
    """
      Special method :
    """
    setattr(self, key, item)

    return
  #
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  #
  def __contains__(self, item):
    """
      Special method :
    """

    if (item in self.__dict__.keys()):
      return True
    else:
      return False
  #
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  #
  def __set_tv_defaults(self):

    cmaps = {}
    cmaps['stki']='viridis'
    cmaps['stkq']='RdGy'
    cmaps['stku']='RdGy'
    cmaps['stkv']='RdGy'
    self.__tv_defaults['cmap']=cmaps

    cscale = {}
    cscale['stki']='linear'
    cscale['stkq']='linear'
    cscale['stku']='linear'
    cscale['stkv']='linear'
    self.__tv_defaults['cscale']=cscale

    factor = {}
    factor['stki']=1.e0
    factor['stkq']=1.e3
    factor['stku']=1.e3
    factor['stkv']=1.e3
    self.__tv_defaults['factor']=factor

    symmetry = {}
    symmetry['stki']=False
    symmetry['stkq']=True
    symmetry['stku']=True
    symmetry['stkv']=True
    self.__tv_defaults['symmetry']=symmetry

    label = {}
    label['stki']=[r'I/I$_{{\rm c}}$', r'I']
    label['stkq']=[r'Q/I$_{{\rm c}}$ [$\times10^{3}$]', r'Q [$\times10^{3}$]']
    label['stku']=[r'U/I$_{{\rm c}}$ [$\times10^{3}$]', r'U [$\times10^{3}$]']
    label['stkv']=[r'V/I$_{{\rm c}}$ [$\times10^{3}$]', r'V [$\times10^{3}$]']
    self.__tv_defaults['label']=label

    self.__tv_defaults['xlabel']=[r'$\lambda$ [px]', r'$\lambda$ [m$\AA$]']

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

  def write(self, fname, vv=3, fmt_type=np.float32):
    self.write_profile(fname, vv=vv, fmt_type=fmt_type)

  def write_profile(self, fname, vv=3, fmt_type=np.float32):

    def set_to_write(obj,v=None):
  
      if (v==None):
        array = np.zeros((4, obj.nw, obj.ny, obj.nx))
        array[0,:,:,:] = np.moveaxis(obj.stki, [0,1,2], [2,0,1]).T * 1.
        array[1,:,:,:] = np.moveaxis(obj.stkq, [0,1,2], [2,0,1]).T * 1.
        array[2,:,:,:] = np.moveaxis(obj.stku, [0,1,2], [2,0,1]).T * 1.
        array[3,:,:,:] = np.moveaxis(obj.stkv, [0,1,2], [2,0,1]).T * 1.
      else:
        array = np.zeros((obj.nx, obj.ny, obj.nw, 4))
        array[:,:,:,0] = np.moveaxis(obj.stki, [0,1,2], [2,0,1]) * 1.
        array[:,:,:,1] = np.moveaxis(obj.stkq, [0,1,2], [2,0,1]) * 1.
        array[:,:,:,2] = np.moveaxis(obj.stku, [0,1,2], [2,0,1]) * 1.
        array[:,:,:,3] = np.moveaxis(obj.stkv, [0,1,2], [2,0,1]) * 1.
  
      return obj.indx, obj.wave, array

    def write_profile3D_v3(fname, data, fmt_type=np.float32):
    
      indx = data[0]
      wave = data[1]
      #CAREFUL, DATA COMES WITH AN ALREADY CHANGED AXIS ORDER!!!!
      # IT MIGHT BE WORTH CHANGING IT?
      # MOREOVER, SO MANY COPIES OF THE SAME THING IS VERY MEMORY CONSUMING!!!
      data4D = data[2]
    
      #_, nw, ny, nx = data4D.shape
      nx, ny, nw, ns = data4D.shape
      #print(data4D.reshape(-1,ns).sum(0))
      data = data4D.flatten()
      del(data4D)
    
      idd = 160904
      v=3
    
      if (fmt_type != np.float32):
        print('Not implemented yet!')
        return np.nan
    
      # We write single precision numbers: 32 bits
      # This is 4 bytes
    
      nnum_to_write = 4. * nw * ny * nx
      nbyts_to_write = nnum_to_write * 4.
    
      #print(nnum_to_write)
    
      # The maximum size of a fortran record is 2Gb
      nbyts_max_frec = 1. * 1024. * 1024. * 1024.
      nnum_max_frec = nbyts_max_frec / 4.
    
      # We set the maximum number of fortran rec to actual_max-10:
      nnum_rec = np.int(nnum_max_frec - 10.)
      nrecs_to_write = np.int32(np.ceil(nnum_to_write/nnum_rec))
    
      # FIRST, WE WRITE A FIRST RECORD WITH THE MANDATORY DATA:
      towrite = np.zeros(9, dtype=fmt_type)
    
      towrite[0] = 3000 + v
      towrite[1] = 3000 - v
      towrite[2] = idd * 1.
      towrite[3] = nrecs_to_write+3.
      towrite[4] = 4.  # NUMBER OF DIMENSIONS
      towrite[5] = nx * 1.
      towrite[6] = ny * 1.
      towrite[7] = nw * 1.
      towrite[8] = 4.  # NUMBER OF STOKES
    
      #OPEN WRITING FILE:
      f=FortranFile(fname,'w')
      # FIRST ACCESS:
      f.write_record(np.float32(towrite))
      # ONGOING ACCESES:
      f.write_record(np.float32(indx))
      f.write_record(np.float32(wave))
      #         SET OF DATA:
      offset = 0
      for it_nnn in range(nrecs_to_write):
        #print('Writing %i of %i' % (it_nnn+1, nrecs_to_write))
        #towrite = data4D[it_nnn,:,:,:].reshape(-1) * 1.
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

    def write_profile3D_v4(fname, data, fmt_type=np.float32):
    
      indx = data[0]
      wave = data[1]
      #CAREFUL, DATA COMES WITH AN ALREADY CHANGED AXIS ORDER!!!!
      # IT MIGHT BE WORTH CHANGING IT?
      # MOREOVER, SO MANY COPIES OF THE SAME THING IS VERY MEMORY CONSUMING!!!
      data4D = data[2]
    
      #_, nw, ny, nx = data4D.shape
      nx, ny, nw, ns = data4D.shape
      #print(data4D.reshape(-1,ns).sum(0))
      data = data4D.flatten()
      del(data4D)
    
      idd = 160904
      v=4
    
      if (fmt_type != np.float32):
        print('Not implemented yet!')
        return np.nan
    
      # We write single precision numbers: 32 bits
      # This is 4 bytes
    
      nnum_to_write = 4. * nw * ny * nx
      nbyts_to_write = nnum_to_write * 4.
    
      #print(nnum_to_write)
    
      # The maximum size of a fortran record is 2Gb
      nbyts_max_frec = 1. * 1024. * 1024. * 1024.
      nnum_max_frec = nbyts_max_frec / 4.
    
      # We set the maximum number of fortran rec to actual_max-10:
      nnum_rec = np.int(nnum_max_frec - 10.)
      nrecs_to_write = np.int32(np.ceil(nnum_to_write/nnum_rec))
    
      # FIRST, WE WRITE A FIRST RECORD WITH THE MANDATORY DATA:
      towrite = np.zeros(9, dtype=fmt_type)
    
      towrite[0] = 3000 + v
      towrite[1] = 3000 - v
      towrite[2] = idd * 1.
      towrite[3] = nrecs_to_write+3.
      towrite[4] = 4.  # NUMBER OF DIMENSIONS
      towrite[5] = nx * 1.
      towrite[6] = ny * 1.
      towrite[7] = nw * 1.
      towrite[8] = 4.  # NUMBER OF STOKES
    
      #OPEN WRITING FILE:
      f=FortranFile(fname,'w')
      # FIRST ACCESS:
      f.write_record(np.float32(towrite))
      # ONGOING ACCESES:
      f.write_record(np.float32(indx))
      f.write_record(np.float32(wave))
      #         SET OF DATA:
      offset = 0
      for it_nnn in range(nrecs_to_write):
        #print('Writing %i of %i' % (it_nnn+1, nrecs_to_write))
        #towrite = data4D[it_nnn,:,:,:].reshape(-1) * 1.
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

    if (vv == 3):
      write_profile3D_v3(fname, set_to_write(self, v=3), fmt_type=fmt_type)
    elif (vv == 4):
      write_profile3D_v4(fname, set_to_write(self, v=4), fmt_type=fmt_type)
    else:
      print('Version %i for profile file not supported!' % i (vv, ))

    return

  def copy(self):

    new_profiles = stk_profile3D(self.nw, self.nx, self.ny)
    new_profiles.set_indx(self.indx)
    new_profiles.set_wave(self.wave)
    new_profiles.set_stki(self.stki)
    new_profiles.set_stkq(self.stkq)
    new_profiles.set_stku(self.stku)
    new_profiles.set_stkv(self.stkv)

    return new_profiles

  def plot(self, **fkwargs):

    plot_profiles([self], **fkwargs)

    return
  #
  def __get_pars(self,pars):

    if (pars[0].lower()=='all'):
      ipars = ['stki', 'stkq', 'stku', 'stkv']
    else:
      ipars = []
      for itp in pars:
        if("stk" in itp):
          to_add = itp
        else:
          to_add = "stk%s" % (itp,)
        if (not(to_add in ipars)):
          ipars.append(to_add)

    return ipars
  #
  # TV
  #
  def tv(self \
      , pars=['all',], waves=[0], fignum=1\
      , axis='p', fkwargs={}, ikwargs=[] \
      , intkwars={}, nearest=False):
    """

      tv method for stk_profile3D class

        inputs:
          [ mandatory ]
            * profiles: one or a list of stk_profile3d class elements.
          [ optional ]
            * pars (default='all'):
            * axis (default='p'):
            * labels (default=None):

            * fkwargs (default={'num':1}): keywords to pass to matplotlib.pyplot subplots call
                                           avoiding ncols, nrows, sharex, squeeze, clear, and sharey
            * pkwargs (default='default matplotlib.pyplot kwargs'): one or a list of dictionaries with
                                           the (matplotlib.pyplot plot standard) keywords to pass to
                                           matplotlib.pyplot plot call.
                                           Length of the list: 1 or the same as the supplied profiles
                                           (mandatory argument) length


    """
    if (type(pars)!=list):
      pars=[pars,]
    if (type(waves)!=list):
      waves=[waves,]

    ipars = self.__get_pars(pars)

    if ( (axis!='p') & (axis!='w') ):
      print("\tAxis must be either 'p' (pixel scale) or 'w' (wavelength scale)")
      return

 #   if (axis=='p'):
 #     xit = np.arange(self.nw, dtype=np.float32)
 #   elif (axis=='w'):
 #     xit = getattr(self, 'wave')

    toshow = np.zeros((len(waves), len(ipars), self.nx, self.ny))

 #   ylabels = []
 #   for itw in len(waves):
 #     if (axis=""):
 #       seedy = ""
 #     else:
 #       seedy = ""
 #   ylabel[itn])


    for itnh, ith in enumerate(waves):
      if ( (ith<0) | ((ith+1)>self.nw)):
        print('\tOut of wavelength range: %i out of [%i-%i]' % (ith, 0, self.nw-1))
        return
      for itnp, itp in enumerate(ipars):
        itpar = getattr(self, itp)
        toshow[itnh, itnp, :, :]=itpar[ith,:,:]*1.

    # Now we start with the plot:
    #
    def show_col(axc, ts, ccmap, factor, label, sym=False, alleq=False):
  
      nts,nx,ny = ts.shape
  
      vmax = np.zeros((nts, ))
      vmin = np.zeros((nts, ))
      for itn in range(nts):
        its = ts[itn,:,:]
        if (np.any(its==its)):
          vmin[itn] = np.percentile(its[its==its], 1) * factor
          vmax[itn] = np.percentile(its[its==its], 99) * factor

      if (alleq==True):
        if (np.any(ts==ts)):
          vmin[:] = np.min(vmin,keepdims=True)
          vmax[:] = np.max(vmax,keepdims=True)
  
      if (sym==True):
        for it_nnn in range(nts):
          vmax[it_nnn] = np.max([np.abs(vmax[it_nnn]), np.abs(vmin[it_nnn])])
          vmin[it_nnn] = -vmax[it_nnn]
  
      for itn in range(nts):
        its = ts[itn,:,:]
        im = axc[itn].imshow(its.T*factor, vmin=vmin[itn], vmax=vmax[itn], cmap=ccmap)
        if (axc[itn].yaxis_inverted()):
          axc[itn].invert_yaxis()
        cbar = pl.colorbar(im, ax=axc[itn], fraction=0.1, shrink=0.9, aspect=30)
        # Following: https://stackoverflow.com/questions/22012096/how-to-set-number-of-ticks-in-plt-colorbar
        tick_locator = ticker.MaxNLocator(nbins=5, min_n_ticks=3)
        cbar.locator = tick_locator
        cbar.update_ticks()
        cbar.set_label(label)
  
      return

    #
    autofigsize=0
    if (not ('figsize' in fkwargs)):
      autofigsize=1
      xs, ys = pl.rcParams['figure.figsize']
      # Get default figure aspect ratio:
      asp_rat_1 = xs / ys
      # Get data aspect ratio:
      asp_rat_2 = (1.*self.nx) / (1.*self.ny)
      # Get toshow aspect ratio:
      asp_rat_3 = (1.*len(ipars)) / (1.*len(waves))
      asp_rat_4 = asp_rat_2 * asp_rat_3
      #
      if (asp_rat_4>1.):
        nxs = xs * asp_rat_4
        nys = ys * 1.
      else:
        nxs = xs * 1.
        nys = ys / asp_rat_4
      nmax = np.max([nxs,nys]) / 10.
      nxs /= nmax
      nys /= nmax
      print(nxs, nys)
      fkwargs['figsize'] = (nxs,nys)
    #
    #
    #
    if (not("num" in fkwargs)):
      fkwargs["num"] = 1

    pl.close(fkwargs["num"])
    fg, ax = pl.subplots(ncols=len(ipars),nrows=len(waves) \
        , squeeze=False, sharex=True, sharey=True, clear=True, **fkwargs)

    #
    #
    #
    itnorm=1
    if (self.stki.mean()<5.):
      itnorm=0

    for itnp, itp in enumerate(ipars):
      show_col(ax[:,itnp], toshow[:,itnp,:,:] \
          , self.__tv_defaults['cmap'][itp] \
          , self.__tv_defaults['factor'][itp] \
          , self.__tv_defaults['label'][itp][itnorm] \
          , sym=self.__tv_defaults['symmetry'][itp])
      # Remove some ticklabels:
      # Y:
      for itn in range(len(waves)):
        ax[itn,itnp].set_ylabel('y [px]')#\n%s' ylabel[itn])
      # X:
      ax[-1,itnp].set_xlabel('x [px]')

#    if (autofigsize==1):
#      pl.rcParams['figure.figsize']=(xs,ys)

    fg.tight_layout()
    pl.show()

    return

    #
    # TV end.
    #

  def trim_profiles(self, trim_dims):

    if (len(trim_dims)!=3):
      print(' ')
      print(' You must supply values for the three dimensions')
      print(' ')
      return

    new_dims = []
    new_tdims = []
    selfs = [self.nw, self.nx, self.ny]
    for it in range(3):
      if (len(trim_dims[it])!=2):
        print(' ')
        print(' You must supply the first and last pixel for each direction')
        print(' ')
        return
      it_dim = []

      if (trim_dims[it][0]<0):
        it_dim.append(selfs[it]+trim_dims[it][0])
      else:
        it_dim.append(trim_dims[it][0])
      if (trim_dims[it][1]<0):
        it_dim.append(selfs[it]+trim_dims[it][1]+1)
      else:
        it_dim.append(trim_dims[it][1])

      if (it_dim[0]>it_dim[1]):
        print(' ')
        print(' Lower limit above upper limit!')
        print(' ')
        return
      if ( (it_dim[0]<0) | (it_dim[1]>selfs[it]) ):
        print(' ')
        print(' Out of boundaries!')
        print(' ')
        return
      new_dims.append(it_dim)
      new_tdims.append(it_dim[1]-it_dim[0])

    new_profiles = stk_profile3D(*tuple(new_tdims))
    pars = ['stki', 'stkq', 'stku', 'stkv']
    for itn, itp in enumerate(pars):
      dum = getattr(self, itp)
      ndum = dum[new_dims[0][0]:new_dims[0][1],new_dims[1][0]:new_dims[1][1],new_dims[2][0]:new_dims[2][1]]
      setattr(new_profiles, itp, ndum)
    pars = ['wave', 'indx']
    for itn, itp in enumerate(pars):
      dum = getattr(self, itp)
      ndum = dum[new_dims[0][0]:new_dims[0][1]]
      setattr(new_profiles, itp, ndum)


    return new_profiles

  #
  # stk_profile3D class end.
  #

#
# Compare TV
#
def compare_tv(stokes, pars, waves, fignum=1\
    , fkwargs={}, ikwargs=[] \
    , intkwars={}):

  if (type(stokes)!=list):
    stokes=[stokes,]
  if (type(pars)!=list):
    pars=[pars,]
  if (type(waves)!=list):
    waves=[waves,]

  # Check shape of all the stk_profiles3D:
  arethesame = True
  for it_stk in range(len(stokes)):
    if (it_stk==0):
      ref=stokes[it_stk].shape
    mult=False
    if (ref==stokes[it_stk].shape):
      mult=True
    arethesame = arethesame * mult
  if (arethesame==False):
    print('\tSupplied stk_profiles3D must have the same shape:')
    for it_stk in range(len(stokes)):
      print(stokes[it_stk].shape)
    return

  if (pars[0].lower()=='all'):
    pars = ['stki', 'stkq', 'stku', 'stkv']

  gnx = stokes[0].nx
  gny = stokes[0].ny

  toshow = np.zeros((len(waves), len(pars), gnx*len(stokes), gny))

  for itnp, itp in enumerate(pars):
    for it_stk in range(len(stokes)):
      itpar = getattr(stokes[it_stk], itp)
      for itnh, ith in enumerate(waves):
        if ( (ith<0) | ((ith+1)>stokes[it_stk].nw)):
          print('\tOut of wavelength range: %i out of [%i-%i]' % (ith, 0, stokes[it_stk].nw-1))
          return
        toshow[itnh, itnp, it_stk*gnx:(it_stk+1)*gnx, :]=itpar[ith,:,:]*1.

  # Now we start with the plot:
  #
  def show_col(axc, ts, ccmap, factor, label, sym=False, alleq=False):

    nts,nx,ny = ts.shape

    vmax = np.zeros((nts, ))
    vmin = np.zeros((nts, ))
    for itn in range(nts):
      its = ts[itn,:,:]
      if (np.any(its==its)):
        vmin[itn] = np.percentile(its[its==its], 1) * factor
        vmax[itn] = np.percentile(its[its==its], 99) * factor

    if (alleq==True):
      if (np.any(ts==ts)):
        vmin[:] = np.min(vmin,keepdims=True)
        vmax[:] = np.max(vmax,keepdims=True)

    if (sym==True):
      for it_nnn in range(nts):
        vmax[it_nnn] = np.max([np.abs(vmax[it_nnn]), np.abs(vmin[it_nnn])])
        vmin[it_nnn] = -vmax[it_nnn]

    for itn in range(nts):
      its = ts[itn,:,:]
      im = axc[itn].imshow(its.T*factor, vmin=vmin[itn], vmax=vmax[itn], cmap=ccmap)
      if (axc[itn].yaxis_inverted()):
        axc[itn].invert_yaxis()
      cbar = pl.colorbar(im, ax=axc[itn], fraction=0.1, shrink=0.9, aspect=30)
      # Following: https://stackoverflow.com/questions/22012096/how-to-set-number-of-ticks-in-plt-colorbar
      tick_locator = ticker.MaxNLocator(nbins=5, min_n_ticks=3)
      cbar.locator = tick_locator
      cbar.update_ticks()
      cbar.set_label(label)

    return

  #
  autofigsize=0
  if (not ('figsize' in fkwargs)):
    autofigsize=1
    xs, ys = pl.rcParams['figure.figsize']
    # Get default figure aspect ratio:
    asp_rat_1 = xs / ys
    # Get data aspect ratio:
    asp_rat_2 = (1.*gnx*len(stokes)) / (1.*gny)
    # Get toshow aspect ratio:
    asp_rat_3 = (1.*len(pars)) / (1.*len(waves))
    aps_rat_4 = asp_rat_2 * asp_rat_3
    #
    if (aps_rat_4>1.):
      nxs = xs * aps_rat_4
      nys = ys * 1.
    else:
      nxs = xs * 1.
      nys = ys / aps_rat_4
    print(nxs, nys)
    pl.rcParams.update({'figure.figsize':(nxs,nys)})
  #
  pl.close(fignum)
  fg, ax = pl.subplots(ncols=len(pars),nrows=len(waves), num=fignum \
      , squeeze=False, **fkwargs)

  itnorm=1
  if (stokes[0].stki.mean()<5.):
    itnorm=0

  for itnp, itp in enumerate(pars):
    show_col(ax[:,itnp], toshow[:,itnp,:,:] \
        , stokes[0].__tv_defaults['cmap'][itp] \
        , stokes[0].__tv_defaults['factor'][itp] \
        , stokes[0].__tv_defaults['label'][itp][itnorm] \
        , sym=stokes[0].__tv_defaults['symmetry'][itp])
    # Remove some ticklabels:
    # Y:
    if ( (len(pars)>1) & (itnp>0) ):
      for itn in range(len(waves)):
        ax[itn,itnp].yaxis.set_ticklabels([])
    if (itnp==0):
      for itn in range(len(waves)):
        ax[itn,itnp].set_ylabel(r'y [px]')
    # X:
    if (len(waves)>1):
      for itn in range(len(waves)-1):
        ax[itn,itnp].xaxis.set_ticklabels([])
      ax[-1,itnp].set_xlabel(r'%i$\times$x [px]' % (len(stokes),))
  
  
  fg.tight_layout()
  pl.show()

  if (autofigsize==1):
    pl.rcParams['figure.figsize']=(xs,ys)

  return

  #
  # compare TV end.
  #


def plot_profiles(profiles \
    , pars=['all'], axis='p', labels=[] \
    , rangex=[], rangey=[] \
    , itx=[0,], ity=[0,] \
    , fkwargs={} \
    , pkwargs={}):
  """

    plot_profiles method:

      inputs:
        [ mandatory ]
          * profiles: one or a list of stk_profile3d class elements.
        [ optional ]
          * pars (default='all'):
          * axis (default='p'):
          * labels (default=None):

          * rangex (default='as given from profiles'): list of two elements giving the minimum
                                         and maximum limits for the x axis

          * fkwargs (default={'num':1}): keywords to pass to matplotlib.pyplot subplots call
                                         avoiding ncols, nrows, sharex, squeeze, clear, and sharey
          * pkwargs (default='default matplotlib.pyplot kwargs'): one or a list of dictionaries with
                                         the (matplotlib.pyplot plot standard) keywords to pass to
                                         matplotlib.pyplot plot call.
                                         Length of the list: 1 or the same as the supplied profiles
                                         (mandatory argument) length

  """

  if (type(profiles)!=list):
    profiles=[profiles,]
  #
  if (type(itx)!=list):
    itx=[itx,]
  if (type(ity)!=list):
    ity=[ity,]
  if (type(labels)!=list):
    labels=[labels,]
  if (type(pars)!=list):
    pars=[pars,]
  if (type(pkwargs) != list):
    pkwargs=[pkwargs,]

  if ( (len(rangex)!=0) & (len(rangex)!=2) ):
    print('rangex must be a list of two elements or not supplied')
    return

  show_labels = True
  if (len(labels)==0):
    show_labels = False
    labels = [''] * len(profiles)

  ipars = profiles[0].__get_pars(pars)

  if ( (len(labels)!=len(profiles)) ):
    print("Length of labels must be equal to profiles' length")
    print("\tLabels length: %i" % (len(labels),))
    print("\tLabels profiles: %i" % (len(profiles),))
    return

  if (len(pkwargs)==1):
    pkwargs = pkwargs * len(profiles)
  #if (len(skwargs)==1):
  #  skwargs = skwargs * len(profiles)

  #
  # Plot itself:
  #
  #
  # Number of columns and rows:
  #
  sqrtabove = np.int32(np.ceil(np.sqrt(len(ipars))))
  axestorem = 0
#  if (len(pars)%2 != 0):
  ncols=np.max([1,(sqrtabove//2)*2])
  nrows=np.int32(np.ceil((1.*len(ipars))/(1.*ncols)))
  axestorem = ncols * nrows - len(ipars)
#  else:
#    ncols=sqrtabove
#    nrows=np.int32(np.ceil(len(pars)/np.float(ncols)))

  if (not("num" in fkwargs)):
    fkwargs["num"] = 1

  pl.close(fkwargs["num"])
  fg,ax=pl.subplots(ncols=ncols,nrows=nrows,sharex=True\
      ,squeeze=False, clear=True, **fkwargs)
  fax = ax.flatten()

  lrangex=1.e99
  urangex=-1.e99
  pltcnt = -1
  for itn, self in enumerate(profiles):

    if (axis=='p'):
      xtoplot = np.arange(self.nw)
    elif(axis=='w'):
      xtoplot = self.wave * 1.
#    else:
#      print("axis='p' or axis='w'")
#      return

    if (np.min(xtoplot)<lrangex):
      lrangex = np.min(xtoplot)
    if (np.max(xtoplot)>urangex):
      urangex = np.max(xtoplot)
  
    for it_nnn in range(len(itx)):

      pltcnt += 1

      for cnt, itp in enumerate(ipars):

        itdata = getattr(self, itp)
        if (show_labels):
          fax[cnt].plot(xtoplot \
              , itdata[:,itx[it_nnn],ity[it_nnn]], label=labels[itn] \
              , **pkwargs[itn])
        else:
          fax[cnt].plot(xtoplot \
              , itdata[:,itx[it_nnn],ity[it_nnn]], **pkwargs[itn] \
              )
        del(itdata)

  #
  #
  #
  ylabel_dict = {}
  #
  if (np.mean(profiles[0].stki)<5.):
    ylabel_dict['stki'] = r'I/I$_{{\rm c}}$'
    ylabel_dict['stkq'] = r'Q/I$_{{\rm c}}$'
    ylabel_dict['stku'] = r'U/I$_{{\rm c}}$'
    ylabel_dict['stkv'] = r'V/I$_{{\rm c}}$'
  else:
    ylabel_dict['stki'] = r'I$'
    ylabel_dict['stkq'] = r'Q$'
    ylabel_dict['stku'] = r'U$'
    ylabel_dict['stkv'] = r'V$'
  for cnt, itp in enumerate(ipars):
    fax[cnt].set_ylabel(ylabel_dict[itp])
  #
  #
  #
  for itn in range(ncols-axestorem):
    if (axis == 'p'):
      ax[-1, itn].set_xlabel(r'$\lambda$ [px]')
    elif (axis == 'w'):
      ax[-1, itn].set_xlabel(r'$\lambda$ [m$\AA$]')

  for itn1 in range(ncols):
    for itn2 in range(nrows-1):
      if (itn2==nrows-2):
        if (axestorem>0):
          if (itn1+axestorem>=ncols):
            pl.delaxes(ax[itn2+1, itn1])
            if (axis == 'p'):
              ax[itn2, itn1].set_xlabel(r'$\lambda$ [px]')
            elif (axis == 'w'):
              ax[itn2, itn1].set_xlabel(r'$\lambda$ [m$\AA$]')
  #
  #
  #
  if (len(rangex)==2):
    lrangex=rangex[0]*1.
    urangex=rangex[1]*1.

    for it_xxx in range(3):
      for it_yyy in range(3):
        ax[it_xxx, it_yyy].set_xlim(lrangex, urangex)
  #
  #
  #
  fg.tight_layout()
  if (show_labels):
    ax[-1,-1].legend(shadow=True, fancybox=True \
        , ncol=np.int(np.round(np.sqrt(len(profiles)))))
  pl.draw()
#, bbox_to_anchor=(0.5,-1.) \

  return





################################################################################
# Model:
################################################################################

def read_model(fname, fmt_type=np.float32, devel=False):

##
## Develpoing:
  def dev_read_model3D_v3(ofile, nrec, dims, fmt_type, verbose=False):
  
    npar, nx, ny, nz = np.int64(dims)
  
    if (np.abs(npar-13.) > 1.e-3):
      print('Something is wrong with the file %s' % fname)
      print('\tIs it a 3D model file?')
      return np.nan
  
    model3D = atm_model3D(nx, ny, nz)

    ipar = 0
    loop = False
    soff = 0
    roff = 0

    ntot = np.int((1.*nx)*ny*nz)
    while (ipar<npar):

      if (loop==True):
        roff = rlast_acc
        soff = 0
      else:
        soff = soff%ntot
        roff=0

      to_store = np.zeros(ntot, dtype=np.float32)
      while (soff<ntot):
        #print("\t", ipar, soff, ntot)

        if (loop==False):
          tmp = ofile.read_record(dtype=fmt_type)

        rlast_acc = np.min([tmp.size, ntot-soff])
        ssize = rlast_acc - roff
        last_acc = soff + ssize

        loop=False
        if (rlast_acc!=tmp.size):
          loop=True

        #print("\t", soff, last_acc, ntot)
        #print("\t\t", roff, rlast_acc, tmp.size)
        #stop()

        to_store[soff:last_acc] = tmp[roff:rlast_acc]*1.

        soff = last_acc * 1
        roff = rlast_acc % tmp.size

      #print(soff, ntot, loop, last_acc, ipar)

      #stop()

      if (ipar==0):
        model3D.set_tem(to_store.reshape(nx,ny,nz))
      elif (ipar==1):
        model3D.set_pg(to_store.reshape(nx,ny,nz))
      elif (ipar==2):
        model3D.set_rho(to_store.reshape(nx,ny,nz))
      elif (ipar==3):
        model3D.set_bx(to_store.reshape(nx,ny,nz))
      elif (ipar==4):
        model3D.set_by(to_store.reshape(nx,ny,nz))
      elif (ipar==5):
        model3D.set_bz(to_store.reshape(nx,ny,nz))
      elif (ipar==6):
        model3D.set_vz(to_store.reshape(nx,ny,nz))
      elif (ipar==9):
        model3D.set_tau(to_store.reshape(nx,ny,nz))
      elif ((model3D.full==True) & (ipar==7)):
        model3D.set_pel(to_store.reshape(nx,ny,nz))
      elif ((model3D.full==True) & (ipar==8)):
        model3D.set_mw(to_store.reshape(nx,ny,nz))
      elif ((model3D.full==True) & (ipar==10)):
        model3D.set_x(to_store.reshape(nx,ny,nz))
      elif ((model3D.full==True) & (ipar==11)):
          model3D.set_y(to_store.reshape(nx,ny,nz))
      elif (ipar==12):
        model3D.set_z(to_store.reshape(nx,ny,nz))

      # Update physical parameter:
      ipar +=1

    del(to_store)
    if (verbose==True):
      print('reading model: %s' % (fname, ))
  
    return model3D
##
## Develpoing.

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
    if (devel==False):
      model3D = read_model3D_v3(f,nrec,first_rec[5:],fmt_type)
    else:
      model3D = dev_read_model3D_v3(f,nrec,first_rec[5:],fmt_type)
  else:
    print('Version %i for model atmosphere file not supported!' % np.int(vers))
    return np.nan

  f.close()

  return model3D

class atm_model3D(object):

  def __init__(self, nx, ny, nz, full=True):

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
    self.tau = np.zeros((self.nx, self.ny, self.nz))
    self.z = np.zeros((self.nx, self.ny, self.nz))
    if (full==True):
      self.pel = np.zeros((self.nx, self.ny, self.nz))
      self.mw = np.zeros((self.nx, self.ny, self.nz))
      self.x = np.zeros((self.nx, self.ny, self.nz))
      self.y = np.zeros((self.nx, self.ny, self.nz))
    else:
      self.pel = np.zeros((0,0,0,))
      self.mw = np.zeros((0,0,0,))
      self.x = np.zeros((0,0,0,))
      self.y = np.zeros((0,0,0,))
    self.shape = self.tem.shape
    self.full = full

    # For tv:
    self.__tv_defaults = {}
    self.__set_tv_defaults()

    return
  #
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  #
  def __repr__(self):
    """
      Special method 'repr':
    """

    return "\n\tFirtez-dz Model atmosphere class."
  #
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  #
  def __add__(self, atmcls_to_add):
    """
      Special method :
    """

    #
    # Check alignment:
    if ( (self.nx!=atmcls_to_add.nx) | (self.ny!=atmcls_to_add.ny) ):
      print("")
      print("\tError!")
      print("\tModel parameters to add are not aligned!")
      print("")
      return np.nan

    out = atm_model3D(self.nx, self.ny, self.nz+stkcls_to_add.nz)


    pars = ["tem", "pg", "rho", "bx", "by", "bz", "vz", "tau", "z", "pel", "mw", "x", "y"]

    for itp in pars:
      if (not itp in out):
        #
        # Second 
        from pdb import set_trace as stop
        stop()
      out[itp] = np.concatenate([self[itp],atmcls_to_add[itp]], axis=2)
        

    return out
  #
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  #
  def __iadd__(self, atmcls_to_add):
    """
      Special method :
    """

    #
    # Check alignment:
    if ( (self.nx!=atmcls_to_add.nx) | (self.ny!=atmcls_to_add.ny) ):
      print("")
      print("\tError!")
      print("\tModel parameters to add are not aligned!")
      print("")
      return np.nan

    pars = ["tem", "pg", "rho", "bx", "by", "bz", "vz", "tau", "z", "pel", "mw", "x", "y"]

    for itp in pars:
      if (not itp in self):
        #
        # Second 
        from pdb import set_trace as stop
        stop()
      self[itp] = np.concatenate([self[itp],atmcls_to_add[itp]], axis=2)
 
    self.shape = self.tem.shape

    return self
  #
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  #
  def __getitem__(self, key):
    """
      Special method :
    """
    return getattr(self, key)
  #
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  #
  def __setitem__(self, key, item):
    """
      Special method :
    """
    setattr(self, key, item)

    return
  #
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  #
  def __contains__(self, item):
    """
      Special method :
    """

    if (item in self.__dict__.keys()):
      return True
    else:
      return False
  #
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  #
  def __set_tv_defaults(self):

    cmaps = {}
    cmaps['tem']='gist_heat'
    cmaps['pg']='viridis'
    cmaps['rho']='viridis'
    cmaps['bx']='RdGy'
    cmaps['by']='RdGy'
    cmaps['bz']='RdGy'
    cmaps['vz']='seismic'
    cmaps['pel']='viridis'
    cmaps['mw']='viridis'
    cmaps['tau']='viridis'
    cmaps['x']='viridis'
    cmaps['y']='viridis'
    cmaps['z']='viridis'
    self.__tv_defaults['cmap']=cmaps

    cscale = {}
    cscale['tem']='linear'
    cscale['pg']='log'
    cscale['rho']='log'
    cscale['bx']='linear'
    cscale['by']='linear'
    cscale['bz']='linear'
    cscale['vz']='linear'
    cscale['pel']='log'
    cscale['mw']='linear'
    cscale['tau']='linear'
    cscale['x']='linear'
    cscale['y']='linear'
    cscale['z']='linear'
    self.__tv_defaults['cscale']=cscale

    factor = {}
    factor['tem']=1.e-3
    factor['pg']=1.e0
    factor['rho']=1.e0
    factor['bx']=1.e-3
    factor['by']=1.e-3
    factor['bz']=1.e-3
    factor['vz']=1.e-5
    factor['pel']=1.e0
    factor['mw']=1.e0
    factor['tau']=1.e0
    factor['x']=1.e-3
    factor['y']=1.e-3
    factor['z']=1.e-3
    self.__tv_defaults['factor']=factor

    symmetry = {}
    symmetry['tem']=False
    symmetry['pg']=False
    symmetry['rho']=False
    symmetry['bx']=True
    symmetry['by']=True
    symmetry['bz']=True
    symmetry['vz']=True
    symmetry['pel']=False
    symmetry['mw']=False
    symmetry['tau']=False
    symmetry['x']=False
    symmetry['y']=False
    symmetry['z']=False
    self.__tv_defaults['symmetry']=symmetry

    label = {}
    label['tem']=r'T [Kk]'
    label['pg']=r'P$_{g}$ [dyn/cm$^2$]'
    label['rho']=r'$\rho$ [gr/$cm^3$]'
    label['bx']=r'B$_{x}$ [KG]'
    label['by']=r'B$_{y}$ [KG]'
    label['bz']=r'B$_{z}$ [KG]'
    label['vz']=r'v$_{z}$ [km/s]'
    label['pel']=r'P$_{el}$ [dyn?/cm$^2$]'
    label['mw']=r'$\mu$'
    label['tau']=r'$\lg\tau_{5}$'
    label['x']=r'x [Mm]'
    label['y']=r'y [Mm]'
    label['z']=r'z [Mm]'
    self.__tv_defaults['label']=label

    return

  def set_attribute(self,narr,itpar,axis):
    targetarr = getattr(self, itpar)
    if (axis==-1):
      if (np.shape(narr)==np.shape(targetarr)):
        narrdim = targetarr * 0. + narr * 1.
      else:
        print('\tSupplied array must matched target array:')
        print('\t Supplied shape:', np.shape(narr))
        print('\t Target shape:', np.shape(targetarr))
        return
    else:
      if (axis==0):
        if (np.shape(targetarr[:,0,0])==np.shape(narr)):
          narrdim = narr[:,None,None] * np.ones((self.ny \
              ,self.nz,))[None,:,:]
      elif (axis==1):
        if (np.shape(targetarr[0,:,0])==np.shape(narr)):
          narrdim = narr[None,:,None] * np.ones((self.nx \
              ,self.nz,))[:,None,:]
      elif (axis==2):
        if (np.shape(targetarr[0,0,:])==np.shape(narr)):
          narrdim = narr[None,None,:] * np.ones((self.nx \
              ,self.ny,))[:,:,None]
      else:
        print('\tAxis must be in between 0 and 2 (both inclusive)')
        return
    setattr(self, itpar, narrdim)
    return

  def set_tem(self,array,axis=-1):
    self.set_attribute(array, 'tem', axis)
    return

  def set_pg(self,array, axis=-1):
    self.set_attribute(array, 'pg', axis)
    return

  def set_rho(self,array, axis=-1):
    self.set_attribute(array, 'rho', axis)
    return

  def set_bx(self,array, axis=-1):
    self.set_attribute(array, 'bx', axis)
    return

  def set_by(self,array, axis=-1):
    self.set_attribute(array, 'by', axis)
    return

  def set_bz(self,array, axis=-1):
    self.set_attribute(array, 'bz', axis)
    return

  def set_vz(self,array, axis=-1):
    self.set_attribute(array, 'vz', axis)
    return

  def set_pel(self,array, axis=-1):
    self.set_attribute(array, 'pel', axis)
    return

  def set_mw(self,array, axis=-1):
    self.set_attribute(array, 'mw', axis)
    return

  def set_tau(self,array, axis=-1):
    self.set_attribute(array, 'tau', axis)
    return

  def set_x(self,array, axis=-1):
    self.set_attribute(array, 'x', axis)
    return

  def set_y(self,array, axis=-1):
    self.set_attribute(array, 'y', axis)
    return

  def set_z(self,array, axis=-1):
    self.set_attribute(array, 'z', axis)
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
    self.set_tau(array[9,:,:,:])
    if (self.full==True):
      self.set_x(array[10,:,:,:])
      self.set_y(array[11,:,:,:])
      self.set_pel(array[7,:,:,:])
      self.set_mw(array[8,:,:,:])
    self.set_z(array[12,:,:,:])

    return

  def write(self, fname, vv=3, fmt_type=np.float32, verbose=False):
    self.write_model(fname, vv=vv, fmt_type=fmt_type, verbose=verbose)
  def write_model(self, fname, vv=3, fmt_type=np.float32, verbose=False):

    def set_to_write(obj,v=None):
  
      if (v==None):
        array = np.zeros((13, obj.nz, obj.ny, obj.nx))
  
        array[0,:,:,:] = obj.tem.T * 1.
        array[1,:,:,:] = obj.pg.T * 1.
        array[2,:,:,:] = obj.rho.T * 1.
        array[3,:,:,:] = obj.bx.T * 1.
        array[4,:,:,:] = obj.by.T * 1.
        array[5,:,:,:] = obj.bz.T * 1.
        array[6,:,:,:] = obj.vz.T * 1.
        array[9,:,:,:] = obj.tau.T * 1.
        array[12,:,:,:] = obj.z.T * 1.
        if (obj.full==True):
          array[7,:,:,:] = obj.pel.T * 1.
          array[8,:,:,:] = obj.mw.T * 1.
          array[10,:,:,:] = obj.x.T * 1.
          array[11,:,:,:] = obj.y.T * 1.
      elif(v==3):
        array = np.zeros((13, obj.nx, obj.ny, obj.nz))
  
        array[0,:,:,:] = obj.tem * 1.
        array[1,:,:,:] = obj.pg * 1.
        array[2,:,:,:] = obj.rho * 1.
        array[3,:,:,:] = obj.bx * 1.
        array[4,:,:,:] = obj.by * 1.
        array[5,:,:,:] = obj.bz * 1.
        array[6,:,:,:] = obj.vz * 1.
        array[9,:,:,:] = obj.tau * 1.
        array[12,:,:,:] = obj.z * 1.
        if (obj.full==True):
          array[7,:,:,:] = obj.pel * 1.
          array[8,:,:,:] = obj.mw * 1.
          array[10,:,:,:] = obj.x * 1.
          array[11,:,:,:] = obj.y * 1.
  
      return array

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
      nbyts_max_frec = 1. * 1024. * 1024. * 1024.
      nnum_max_frec = nbyts_max_frec / 4.
    
      # We set the maximum number of fortran rec to actual_max-10:
      nnum_rec = np.int(nnum_max_frec - 10.)
      nrecs_to_write = np.int32(np.ceil(nnum_to_write/nnum_rec))
    
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

    if (verbose==True):
      print('Writing model: %s' % (fname,))
    if (vv == 3):
      write_model3D_v3(fname, set_to_write(self, v=vv), fmt_type=fmt_type)
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

  def copy(self):

    new_model = atm_model3D(self.nx, self.ny, self.nz)
    new_model.set_tem(self.tem)
    new_model.set_pg(self.pg)
    new_model.set_rho(self.rho)
    new_model.set_bx(self.bx)
    new_model.set_by(self.by)
    new_model.set_bz(self.bz)
    new_model.set_vz(self.vz)
    new_model.set_pel(self.pel)
    new_model.set_mw(self.mw)
    new_model.set_tau(self.tau)
    new_model.set_x(self.x)
    new_model.set_y(self.y)
    new_model.set_z(self.z)

    return new_model

  def expand_model(self, new_dims):

    if ( len(new_dims) != 3 ):
      print('')
      print(' New dimensions must be a list of length 3!...')
      print(' ... not: ', new_dims)
      print('')
      return

    new_model = atm_model3D(*tuple(new_dims))
    dx = new_model.nx-self.nx
    dy = new_model.ny-self.ny
    dz = new_model.nz-self.nz
    pars = ['tem', 'pg', 'rho', 'bx', 'by', 'bz', 'vz' \
        , 'pel', 'mw', 'tau', 'x', 'y', 'z']
    for itn, itp in enumerate(pars):
      #
      dum = getattr(self, itp)
      #
      amidone = False
      if ( (self.nx != 1) & (itp=='x') ):
        ndum = np.pad(dum, ((0,0), (0,dy), (0,dz) ), mode='edge')
        sdx = ndum[-1,:,:]-ndum[-2,:,:]
        sdx = sdx[None,:,:] * (1.+np.arange(dx))[:,None,None] \
            +ndum[-1,:,:][None,:,:]
        ndum = np.concatenate([ndum, sdx], axis=0)
        amidone = True
      if ( (self.ny != 1) & (itp=='y') ):
        ndum = np.pad(dum, ((0,dx), (0,0), (0,dz) ), mode='edge')
        sdy = ndum[:,-1,:]-ndum[:,-2,:]
        sdy = sdy[:,None,:] * (1.+np.arange(dy))[None,:,None] \
            +ndum[:,-1,:][:,None,:]
        ndum = np.concatenate([ndum, sdy], axis=1)
        amidone = True
      if ( (self.nz != 1) & (itp=='z') ):
        ndum = np.pad(dum, ((0,dx), (0,dy), (0,0) ), mode='edge')
        sdz = ndum[:,:,-1]-ndum[:,:,-2]
        sdz = sdz[:,:,None] * (1.+np.arange(dz))[None,None,:] \
            +ndum[:,:,-1][:,:,None]
        ndum = np.concatenate([ndum, sdz], axis=2)
        amidone = True
      if (amidone == False):
        ndum = np.pad(dum, ((0,dx), (0,dy), (0,dz) ), mode='edge')
      #
      setattr(new_model, itp, ndum)
      #

    return new_model

  def trim_model(self, trim_dims):

    if (len(trim_dims)!=3):
      print(' ')
      print(' You must supply values for the three dimensions')
      print(' ')
      return

    new_dims = []
    new_tdims = []
    selfs = [self.nx, self.ny, self.nz]
    for it in range(3):
      if (len(trim_dims[it])!=2):
        print(' ')
        print(' You must supply the first and last pixel for each direction')
        print(' ')
        return
      it_dim = []

      if (trim_dims[it][0]<0):
        it_dim.append(selfs[it]+trim_dims[it][0])
      else:
        it_dim.append(trim_dims[it][0])
      if (trim_dims[it][1]<0):
        it_dim.append(selfs[it]+trim_dims[it][1]+1)
      else:
        it_dim.append(trim_dims[it][1])

      if (it_dim[0]>it_dim[1]):
        print(' ')
        print(' Lower limit above upper limit!')
        print(' ')
        return
      if ( (it_dim[0]<0) | (it_dim[1]>selfs[it]) ):
        print(' ')
        print(' Out of boundaries!')
        print(' ')
        return
      new_dims.append(it_dim)
      new_tdims.append(it_dim[1]-it_dim[0])
#      print(new_dims)
#      print(new_tdims)

    new_model = atm_model3D(*tuple(new_tdims))
    pars = ['tem', 'pg', 'rho', 'bx', 'by', 'bz', 'vz' \
        , 'pel', 'mw', 'tau', 'x', 'y', 'z']
    for itn, itp in enumerate(pars):
      dum = getattr(self, itp)
      ndum = dum[new_dims[0][0]:new_dims[0][1],new_dims[1][0]:new_dims[1][1],new_dims[2][0]:new_dims[2][1]]
      setattr(new_model, itp, ndum)

    return new_model



  def plot(self, **fkwargs):

    plot_models([self], **fkwargs)

    return

  #
  # TV
  #
  def tv(self, pars, heights, axis='z', fignum=1\
      , fkwargs={}, ikwargs=[] \
      , intkwars={}, nearest=False):

    if (type(pars)!=list):
      pars=[pars,]
    if (type(heights)!=list):
      heights=[heights,]

    if ( (axis!='z') & (axis!='t') ):
      print("\tAxis must be either 'z' (z scale) or 't' (tau scale)")
      return

    if (axis=='z'):
      xit = getattr(self, 'z')
    elif (axis=='t'):
      xit = getattr(self, 'tau')

    ylabels = []
    for itn in range(len(heights)):
      if (axis=="z"):
        prefixy = r"z="
        suffixy = r"[Mm]"
        factory = 1.e-3
      else:
        prefixy = r"$\log\tau_{5}=$"
        suffixy = r""
        factory = 1.e0
      ylabels.append("%s %.2f %s" % (prefixy,heights[itn]*factory,suffixy,))

    toshow = np.zeros((len(heights), len(pars), self.nx, self.ny))

    for itnh, ith in enumerate(heights):

      posclose = np.argmin(np.abs(xit-ith), axis=2)

      if (nearest==False):
        posmin = posclose - 1
        posmax = posclose + 2
  
        for itnp, itp in enumerate(pars):
          itpar = getattr(self, itp)
          for itx in range(self.nx):
            for ity in range(self.ny):
              amin = np.max([0,posmin[itx,ity]])
              amax = np.min([self.nz,posmax[itx,ity]])
              if (amin==0):
                amax = np.max([amax, 3])
              if (amax==self.nz):
                amin = np.min([amin, self.nz-3])
              fint = interp1d(xit[itx,ity,amin:amax] \
                  , itpar[itx,ity,amin:amax] \
                  , bounds_error=False, fill_value=np.nan \
                  , **intkwars)
              toshow[itnh, itnp, itx, ity] = fint(ith)
      else:
        for itnp, itp in enumerate(pars):
          itpar = getattr(self, itp)
          for itx in range(self.nx):
            for ity in range(self.ny):
              toshow[itnh, itnp, itx, ity] = itpar[itx,ity,posclose[itx,ity]]

    # Now we start with the plot:
    #
    def show_col(axc, ts, ccmap, factor, label, sym=False, alleq=False):
  
      nts,nx,ny = ts.shape
  
      vmax = np.zeros((nts, ))
      vmin = np.zeros((nts, ))
      for itn in range(nts):
        its = ts[itn,:,:]
        if (np.any(its==its)):
          vmin[itn] = np.percentile(its[its==its], 1) * factor
          vmax[itn] = np.percentile(its[its==its], 99) * factor

      if (alleq==True):
        if (np.any(ts==ts)):
          vmin[:] = np.min(vmin,keepdims=True)
          vmax[:] = np.max(vmax,keepdims=True)
  
      if (sym==True):
        for it_nnn in range(nts):
          vmax[it_nnn] = np.max([np.abs(vmax[it_nnn]), np.abs(vmin[it_nnn])])
          vmin[it_nnn] = -vmax[it_nnn]
  
      for itn in range(nts):
        its = ts[itn,:,:]
        im = axc[itn].imshow(its.T*factor, vmin=vmin[itn], vmax=vmax[itn], cmap=ccmap)
        if (axc[itn].yaxis_inverted()):
          axc[itn].invert_yaxis()
        cbar = pl.colorbar(im, ax=axc[itn], fraction=0.1, shrink=0.9, aspect=30)
        tick_locator = ticker.MaxNLocator(nbins=5, min_n_ticks=3)
        cbar.locator = tick_locator
        cbar.update_ticks()
        cbar.set_label(label)
  
      return

    #
    #
    autofigsize=0
    if (not ('figsize' in fkwargs)):
      autofigsize=1
      xs, ys = pl.rcParams['figure.figsize']
      # Get default figure aspect ratio:
      asp_rat_1 = xs / ys
      # Get data aspect ratio:
      asp_rat_2 = (1.*self.nx) / (1.*self.ny)
      # Get toshow aspect ratio:
      asp_rat_3 = (1.*len(pars)) / (1.*len(heights))
      aps_rat_4 = asp_rat_2 * asp_rat_3
      #
      if (aps_rat_4>1.):
        nxs = xs * aps_rat_4
        nys = ys * 1.
      else:
        nxs = xs * 1.
        nys = ys / aps_rat_4
      print(nxs, nys)
      pl.rcParams.update({'figure.figsize':(nxs,nys)})
    #

    pl.close(fignum)
    fg, ax = pl.subplots(ncols=len(pars),nrows=len(heights), num=fignum \
        , squeeze=False, sharex=True, sharey=True, **fkwargs)
    for itnp, itp in enumerate(pars):
      show_col(ax[:,itnp], toshow[:,itnp,:,:] \
          , self.__tv_defaults['cmap'][itp] \
          , self.__tv_defaults['factor'][itp] \
          , self.__tv_defaults['label'][itp] \
          , sym=self.__tv_defaults['symmetry'][itp])
      # Remove some ticklabels:
      # Y:
      if (itnp==0):
        for itn in range(len(heights)):
          ax[itn,itnp].set_ylabel('%s\ny [px]' % (ylabels[itn],))
      # X:
      ax[-1,itnp].set_xlabel('x [px]')
    
    fg.tight_layout()
    pl.show()

    if (autofigsize==1):
      pl.rcParams['figure.figsize']=(xs,ys)

    return
  #
  # END atm_model3D class 
  #

#
# START: plot models:
#
def plot_models(models, pars=['all'],fnum=1,itx=[0,],ity=[0,], axis='t', labels=[] \
    , rangex=[], zorigin=False, displaytau=True \
    , skwargs={}, pargs=(), pkargs={}, fkwargs={}):

  # Define some functions for this particular function:

  def get_model_lims(xx, lx, ux, yy, current_lims):

    new_lims = current_lims * 1.

    ww = np.where( (xx-lx) < 0.)[0]
    it_pxl = ww[np.argmax(xx-lx)]
    ww = np.where( (xx-ux) > 0.)[0]
    it_pxu = ww[np.argmin(xx[ww]-ux)]
          #
    low=np.min([it_pxl, it_pxu])
    upp=np.max([it_pxl, it_pxu])
    test_yl = np.min(yy[low:upp+1])
    test_ym = np.mean(np.abs(yy[low:upp+1]))
    test_yu = np.max(yy[low:upp+1])

    if (test_yl<0.):
      test_yl=test_yl - 0.2 * test_ym#*1.1
    if (test_yu<0.):
      test_yu=test_yu + 0.2 * test_ym#*0.9

    if (test_yl>0.):
      test_yl=test_yl - 0.2 * test_ym#*0.9
    if (test_yu>0.):
      test_yu=test_yu + 0.2 * test_ym#*1.1

    if (test_yl < new_lims[0]):
      new_lims[0] = test_yl * 1.
    if (test_yu > new_lims[1]):
      new_lims[1] = test_yu * 1.

    return new_lims

  def oplot_tau_axis(iax, taus, zs, tauvalues, yval, color):
  
  
    xlims = iax.get_xlim()

    ftauz = interp1d(taus, zs, bounds_error=False)#, fill_value="extrapolate"))
    zvalues = ftauz(tauvalues)
  
    axfracs = (np.array(zvalues) - xlims[0]) / (xlims[1] - xlims[0])
  
    iax.plot(axfracs, [yval] * len(axfracs), transform=iax.transAxes, marker='|', color=(1.0,1.,1.), alpha=0.7, linewidth=7, markersize=8)
    iax.plot(axfracs, [yval] * len(axfracs), transform=iax.transAxes, marker='|', color=color)
  
    return

  #
  # Start the function plot_models itself
  #
  # First, some input checks:
  #

  if (type(itx)!=list):
    itx=[itx,]
  if (type(ity)!=list):
    ity=[ity,]
  if (type(models)!=list):
    models=[models,]
  if (type(pars)!=list):
    pars=[pars,]
  if (type(labels)!=list):
    labels=[labels,]
  if (type(skwargs) != list):
    skwargs=[skwargs,]

  if (pars[0]=='all'):
    pars = ['tem', 'pg', 'rho', 'bx', 'by', 'bz', 'vz' \
        , 'pel', 'mw', 'tau', 'x', 'y', 'z']

  if ( (len(rangex)!=0) & (len(rangex)!=2) ):
    print('rangex must be a list of two elements or not supplied')
    return

  show_labels = True
  if (len(labels)==0):
    show_labels = False
    labels = [''] * len(models)

  if ( (len(labels)!=len(models)) ):
    print("Length of labels must be equal to models' length")
    print("\tLabels length: %i" % (len(labels),))
    print("\tLabels models: %i" % (len(models),))
    return

  if (len(skwargs)==1):
    skwargs = skwargs * len(models)

  #
  # Number of columns and rows:
  #
  sqrtabove = np.int16(np.ceil(np.sqrt(len(pars))))
  axestorem = 0
  if (len(pars)%2 != 0):
    ncols=np.max([1,sqrtabove//2*2])
    nrows=np.int16(np.ceil((1.*len(pars))/(1.*ncols)))
    axestorem = ncols * nrows - len(pars)
  else:
    ncols=sqrtabove
    nrows=np.int16(np.ceil(len(pars)/np.float(ncols)))
  #
  # Once input checked, plot atmosphere model:
  #
  #
  autofigsize=0
  if (not ('figsize' in fkwargs)):
    autofigsize=1
    xs, ys = pl.rcParams['figure.figsize']
    nxs = ncols / 1.*np.min([8.,xs])
    nys = np.min([4.,ys]) * nrows
    pl.rcParams.update({'figure.figsize':(nxs,nys)})
  #
  pl.close(fnum)
  fg,ax=pl.subplots(ncols=ncols,nrows=nrows,num=fnum, clear=True, squeeze=False, **fkwargs)
  fax = ax.flatten()
  lrangex=1.e99
  urangex=-1.e99
  zoffsets = np.zeros((len(models),len(itx)))
  used_colors = []

  def plot_par(iax, ix, iy, ifactor, iscale, iskwargs, ilab):
    lin=iax.step(ix, iy * ifactor, label=ilab, **iskwargs)
    iax.set_yscale(iscale)
    return lin[0].get_color()

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

      it_xtoplot = xtoplot[itx[it_nnn],ity[it_nnn],:]

      if ( (axis=='z') & (zorigin==True) ):
        f = interp1d(self.tau[itx[it_nnn],ity[it_nnn],:], it_xtoplot\
            , bounds_error=False, fill_value="extrapolate")
        zoffsets[itn, it_nnn] = - (f(0.)-1.)
        it_xtoplot = it_xtoplot - (f(0.)-1.)

      # We go through the parameters:
      for itpn, itp in enumerate(pars):
        itcolor=plot_par(fax[itpn], it_xtoplot, getattr(self, itp)[itx[it_nnn] \
            ,ity[it_nnn],:], self.__tv_defaults['factor'][itp] \
            , self.__tv_defaults['cscale'][itp], skwargs[itn] \
            , labels[itn])
        if (itpn == 0):
          used_colors.append(itcolor)
        if ( (it_nnn==0) & (itn==0) ):
          fax[itpn].set_ylabel(self.__tv_defaults['label'][itp])

  for itn in range(ncols-axestorem):
    if (axis == 'z'):
      ax[-1, itn].set_xlabel(models[0].__tv_defaults['label']['z'])
    elif (axis == 't'):
      ax[-1, itn].set_xlabel(models[0].__tv_defaults['label']['tau'])
  # Final layout
  fg.tight_layout()

  for itn1 in range(ncols):
    for itn2 in range(nrows-1):
      if (itn2==nrows-2):
        if (axestorem>0):
          if (itn1+axestorem<ncols):
            ax[itn2, itn1].xaxis.set_ticklabels([])
          else:
            pl.delaxes(ax[itn2+1, itn1])
            if (axis == 'z'):
              ax[itn2, itn1].set_xlabel(models[0].__tv_defaults['label']['z'])
            elif (axis == 't'):
              ax[itn2, itn1].set_xlabel(models[0].__tv_defaults['label']['tau'])
      #  else:
      #    ax[itn2, itn1].xaxis.set_ticklabels([])
      #else:
      #  ax[itn2, itn1].xaxis.set_ticklabels([])
          
  if (len(rangex)==2):
    lrangex=rangex[0]*1.
    urangex=rangex[1]*1.

    # store ylims:
    ylims = np.zeros((2,nrows,ncols))
    ylims[0,:,:]=np.inf
    ylims[1,:,:]=-np.inf

    # Set xrange and get the y limits for this range
    cnt = -1
    for it_xxx in range(nrows):
      for it_yyy in range(ncols):
        cnt+=1
        if ((cnt+1)>len(pars)):
          continue
        ax[it_xxx, it_yyy].set_xlim(lrangex, urangex)
        # Ylimits:
        for itn, self in enumerate(models):
          for it_nnn in range(len(itx)):
            #
            if (axis=='z'):
              it_px = self.z[itx[it_nnn],ity[it_nnn],:] * 1.e-3
            elif (axis=='t'):
              it_px = self.tau[itx[it_nnn],ity[it_nnn],:] * 1.
            #
            it_py = getattr(self, pars[cnt])[itx[it_nnn],ity[it_nnn],:] \
                * self.__tv_defaults['factor'][pars[cnt]]

            ylims[:,it_xxx, it_yyy] = get_model_lims(it_px \
                , lrangex, urangex, it_py, ylims[:,it_xxx, it_yyy])

        var = (ylims[1,it_xxx, it_yyy]-ylims[0,it_xxx, it_yyy])/np.max(np.abs(ylims[:,it_xxx, it_yyy]))
        if ( (var > 0.01) & (np.abs(var) != np.inf) & (var==var) ):
          ax[it_xxx, it_yyy].set_ylim(ylims[0,it_xxx, it_yyy], ylims[1,it_xxx, it_yyy])
  # Tau bar and z origin if z plot:
  if (axis=='z'):
    htau = [0.,-1.,-2.,-3.]
    for it_xxx in range(nrows):
      for it_yyy in range(ncols):
        if (zorigin==True):
          ax[it_xxx, it_yyy].axvline(1., color=(0.6,0.6,0.6), linewidth=1., linestyle='--')
    if (displaytau==True):
      cnt=-1
      for itn, self in enumerate(models):
        for it_nnn in range(len(itx)):
          cnt+=1
          it_px = self.z[itx[it_nnn],ity[it_nnn],:] * 1.e-3
          it_py = self.tau[itx[it_nnn],ity[it_nnn],:] * 1.
          #
          oplot_tau_axis(fax[0], it_py, it_px+zoffsets[itn,it_nnn], htau, 0.95-0.05*cnt, used_colors[cnt])
  # Legend:
  if (show_labels):
    if (axestorem>0):
      ax[-2,-1].legend(loc=8,shadow=True, fancybox=True, bbox_to_anchor=(0.5,-1.) \
          , ncol=np.int(np.round(np.sqrt(len(models)))))
    else:
      fax[-1].legend(shadow=True, fancybox=True \
          , ncol=np.int(np.round(np.sqrt(len(models)))))

  if (autofigsize==1):
    pl.rcParams['figure.figsize']=(xs,ys)

  return



################################################################################
# Control file:
################################################################################

#class control(object):
#
#  def __init__(self):
#    self.box = init_box()
#    self.lines = {'hd':'LINES:'}
#    self.mode = {'hd':'MODE:'}
#    self.fileprofile = {'hd':'FILEPROFILE:'}
#    self.filemodel = {'hd':'FILEMODEL:'}
#    self.hydrostatic = {'hd':'HYDROSTATIC:'}
#    self.stokes_setup = {'hd':'STOKES SETUP:'}
#    self.inversion_setup = {'hd':'INVERSION SETUP:'}
#    self.continuum_normalization = {'hd':'CONTINUUM NORMALIZATION:'}
#    self.coupled_inversion = {'hd':'COUPLED INVERSION:'}
#    self.write_respones_function = {'hd':'WRITE RESPONSE FUNCTION:'}
#    self.line_spread_function = {'hd':'LINE SPREAD FUNCTION:'}
#    self.misc_setup = {'hd':'MISC SETUP:'}
#    self.end = {'hd':'END:'}
#  def init_box():
#    rdict = {}
#    rdict['hd']='BOX:'
#    rdict['help']=['nx(i) dx(f)\nny(i) dy(f)\nnz(i) dz(f)']
#    return rdict
#  def init_lines():
#    rdict = {}
#    rdict['hd']='LINES:'
#    rdict['help']=
#    return rdict
#  def init_mode():
#    rdict = {}
#    rdict['hd']='MODE:'
#    rdict['help']=
#    return rdict
#  def init_fileprofile():
#    rdict = {}
#    rdict['hd']='FILEPROFILE:'
#  def init_filemodel():
#    rdict = {}
#    rdict['hd']='FILEMODEL:'
#  def init_hydrostatic():
#    rdict = {}
#    rdict['hd']='HYDROSTATIC:'
#  def init_stokes_setup():
#    rdict = {}
#    rdict['hd']='STOKES SETUP:'
#  def init_inversion_setup():
#    rdict = {}
#    rdict['hd']='INVERSION SETUP:'
#  def init_continuum_normalization():
#    rdict = {}
#    rdict['hd']='CONTINUUM NORMALIZATION:'
#  def init_coupled_inversion():
#    rdict = {}
#    rdict['hd']='COUPLED INVERSION:'
#  def init_write_respones_function():
#    rdict = {}
#    rdict['hd']='WRITE RESPONSE FUNCTION:'
#  def init_line_spread_function():
#    rdict = {}
#    rdict['hd']='LINE SPREAD FUNCTION:'
#  def init_misc_setup():
#    rdict = {}
#    rdict['hd']='MISC SETUP:'

################################################################################
# Response functions:
################################################################################

class stk_rf3D(object):

  def __init__(self, nx, ny, nz, nw, ns, nd):

    self.nx = nx * 1
    self.ny = ny * 1
    self.nz = nz * 1
    self.nw = nw * 1
    self.ns = ns * 1
    self.nd = nd * 1

    dnx = 1
    dny = 1
    dnz = 1
    dnw = 1
    dns = 1
    dnd = 1

    self.z = np.zeros(self.nz)
    self.indx = np.zeros(self.nw)
    self.wave = np.zeros(self.nw)

    # Generate dummy dimension variables, so that they exists but do not...
    # ... take too much memory:
    self.rf_tem = np.zeros((dnx, dny, dnz, dnw, dns))
    self.rf_pg = np.zeros((dnx, dny, dnz, dnw, dns))
    self.rf_rho = np.zeros((dnx, dny, dnz, dnw, dns))
    self.rf_bx = np.zeros((dnx, dny, dnz, dnw, dns))
    self.rf_by = np.zeros((dnx, dny, dnz, dnw, dns))
    self.rf_bz = np.zeros((dnx, dny, dnz, dnw, dns))
    self.rf_vz = np.zeros((dnx, dny, dnz, dnw, dns))
    self.rf_p0 = np.zeros((dnx, dny, dnz, dnw, dns))

    return

  def set_z(self,array):
    self.z = self.z * 0. + array * 1.
    return

  def set_indx(self,array):
    self.indx = self.indx * 0. + array * 1.
    return

  def set_wave(self,array):
    self.wave = self.wave * 0. + array * 1.
    return

  def set_rf_tem(self,array):
    self.rf_tem = np.zeros((self.nx, self.ny, self.nz, self.nw, self.ns))
    self.rf_tem = self.rf_tem * 0. + array * 1.
    return

  def set_rf_pg(self,array):
    self.rf_pg = np.zeros((self.nx, self.ny, self.nz, self.nw, self.ns))
    self.rf_pg = self.rf_pg * 0. + array * 1.
    return

  def set_rf_rho(self,array):
    self.rf_rho = np.zeros((self.nx, self.ny, self.nz, self.nw, self.ns))
    self.rf_rho = self.rf_rho * 0. + array * 1.
    return

  def set_rf_bx(self,array):
    self.rf_bx = np.zeros((self.nx, self.ny, self.nz, self.nw, self.ns))
    self.rf_bx = self.rf_bx * 0. + array * 1.
    return

  def set_rf_by(self,array):
    self.rf_by = np.zeros((self.nx, self.ny, self.nz, self.nw, self.ns))
    self.rf_by = self.rf_by * 0. + array * 1.
    return

  def set_rf_bz(self,array):
    self.rf_bz = np.zeros((self.nx, self.ny, self.nz, self.nw, self.ns))
    self.rf_bz = self.rf_bz * 0. + array * 1.
    return

  def set_rf_vz(self,array):
    self.rf_vz = np.zeros((self.nx, self.ny, self.nz, self.nw, self.ns))
    self.rf_vz = self.rf_vz * 0. + array * 1.
    return

  def set_rf_p0(self,array):
    self.rf_p0 = np.zeros((self.nx, self.ny, self.nz, self.nw, self.ns))
    self.rf_p0 = self.rf_p0 * 0. + array * 1.
    return

  def set_srf(self, array, logic):

    dims = np.array(array.shape)

    if ( (dims.size != 6)\
        | (dims[0] != self.nx)\
        | (dims[1] != self.ny)\
        | (dims[2] != self.nz)\
        | (dims[3] != self.nw)\
        | (dims[4] != self.ns)\
        | (dims[5] != np.int(np.sum(logic)))\
        ):
      print('Wrong format')
      return

    pars = [\
        'tem'
        , 'pg'
        , 'rho'
        , 'bx'
        , 'by'
        , 'bz'
        , 'vz'
        , 'p0']

    cnt = -1
    for itn,itd in enumerate(logic):
      if (np.abs(itd-1.)<1.e-5):
        cnt+=1
        if (itn==0):
          self.set_rf_tem(array[:,:,:,:,:,cnt])
        if (itn==1):
          self.set_rf_pg(array[:,:,:,:,:,cnt])
        if (itn==2):
          self.set_rf_rho(array[:,:,:,:,:,cnt])
        if (itn==3):
          self.set_rf_bx(array[:,:,:,:,:,cnt])
        if (itn==4):
          self.set_rf_by(array[:,:,:,:,:,cnt])
        if (itn==5):
          self.set_rf_bz(array[:,:,:,:,:,cnt])
        if (itn==6):
          self.set_rf_vz(array[:,:,:,:,:,cnt])
        if (itn==7):
          self.set_rf_p0(array[:,:,:,:,:,cnt])

    return

#
def read_rf(fname, itx=-1, fmt_type=np.float32):

  def read_rf3D_v3(ofile, nrec, dims, fmt_type):
  
    nx, ny, nz, nw, ns, nd = np.int64(dims)
  
    # WE ONLY READ A SINGLE X POSITION, OTHERWISE WE CAN...
    # ...RUN OUT OF MEMORY
    rf3D = stk_rf3D(nx, ny, nz, nw, ns, nd)
  
    # IN ORDER TO READ IT, WE KNOW THAT:
    #    1st RECORD: Z
    #    2nd RECORD: INDX
    #    3rd RECORD: WAVE
  
    for it_nnn in range(3):
      tmp = ofile.read_record(dtype=fmt_type)
      if (it_nnn == 0):
        rf3D.set_z(tmp)
      if (it_nnn == 1):
        rf3D.set_indx(tmp)
      if (it_nnn == 2):
        rf3D.set_wave(tmp)
    # logic stores what physical parameters are really used:
    logic = ofile.read_record(dtype=fmt_type)
    npar = np.int64(np.sum(logic))
    #
    #    N   RECORDS: RFS GROUPED BY y, z, w, stokes, pars
    ntotdims = np.int64(nx*ny*nz*nw*ns*npar)

    to_store = np.zeros(ntotdims, dtype=np.float32)
    offset=0
    for it_nnn in range(np.int64(nrec)-5):
      tmp = ofile.read_record(dtype=fmt_type)
      to_store[offset:offset+tmp.size] = tmp*1.
      offset=offset+tmp.size
  
    rf3D.set_srf(to_store.reshape(nx,ny,nz,nw,ns,npar), logic)
    del(to_store)
  
    return rf3D

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

  if ( (np.abs(medv-3000.) > 1.e-3) | (np.abs(fid-18060904.) > 1.e-3) ):
    print('Something is wrong with the file %s' % fname)
    print('\tIs it a 3D response function file?')
    f.close()
    return np.nan

  if (np.abs(vers-3.) < 1.e-3):
    #VERSION 3
    rf3D = read_rf3D_v3(f,nrec,first_rec[5:], fmt_type)
  else:
    print('The version of this RF file is not compatible with this...')
    print('... library version!')
    f.close()
    return np.nan

  f.close()

  return rf3D
#

class logtaus(object):

  def __init__(self, nx, ny, nz, nw):

    self.nx = nx * 1
    self.ny = ny * 1
    self.nz = nz * 1
    self.nw = nw * 1

    self.z = np.zeros(self.nz)
    self.indx = np.zeros(self.nw)
    self.wave = np.zeros(self.nw)

    self.ltau = np.zeros((nx, ny, nz, nw))

    return

  def set_z(self,array):
    self.z = self.z * 0. + array * 1.
    return

  def set_indx(self,array):
    self.indx = self.indx * 0. + array * 1.
    return

  def set_wave(self,array):
    self.wave = self.wave * 0. + array * 1.
    return

  def set_ltau(self,array):
    self.ltau = self.ltau * 0. + array * 1.
    return

def read_logtaus(fname, fmt_type=np.float32):

  def read_ltau_v3(ofile, nrec, dims, fmt_type):
  
    nx, ny, nz, nw = np.int64(dims)
  
    # WE ONLY READ A SINGLE X POSITION, OTHERWISE WE CAN...
    # ...RUN OUT OF MEMORY
    ltaus = logtaus(nx, ny, nz, nw)
  
    # IN ORDER TO READ IT, WE KNOW THAT:
    #    1st RECORD: Z
    #    2nd RECORD: INDX
    #    3rd RECORD: WAVE
  
    for it_nnn in range(3):
      tmp = ofile.read_record(dtype=fmt_type)
      if (it_nnn == 0):
        ltaus.set_z(tmp)
      if (it_nnn == 1):
        ltaus.set_indx(tmp)
      if (it_nnn == 2):
        ltaus.set_wave(tmp)
    #
    #    N   RECORDS: RFS GROUPED BY y, z, w, stokes, pars
    ntotdims = np.int64(nx*ny*nz*nw)

    #to_store = np.zeros(ntotdims, dtype=np.float32)
    ltaus.ltau = ltaus.ltau.reshape(ntotdims)
    offset=0
    for it_nnn in range(np.int64(nrec)-4):
      print(it_nnn)
      tmp = ofile.read_record(dtype=fmt_type)
      ltaus.ltau[offset:offset+tmp.size] = tmp*1.
      offset=offset+tmp.size
  
    #ltaus.set_ltau(to_store.reshape(nx,ny,nz,nw))
    ltaus.ltau = ltaus.ltau.reshape(nx,ny,nz,nw)
    #del(to_store)
  
    return ltaus

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

  if ( (np.abs(medv-3000.) > 1.e-3) | (np.abs(fid-12200904.) > 1.e-3) ):
    print('Something is wrong with the file %s' % fname)
    print('\tIs it a line optical depth file?')
    f.close()
    return np.nan

  if (np.abs(vers-3.) < 1.e-3):
    #VERSION 3
    ltaus = read_ltau_v3(f,nrec,first_rec[5:], fmt_type)
  else:
    print('The version of this RF file is not compatible with this...')
    print('... library version!')
    f.close()
    return np.nan

  f.close()

  return ltaus
#
################################################################################
# Generic reader:
################################################################################

def read(fname, fmt_type=np.float32, devel=False, itx=-1):

  def read_profile3D_v3(ofile, nrec, dims, fmt_type):
  
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

  f = FortranFile(fname, 'r')
  first_rec = f.read_record(dtype=fmt_type)
  f.close()

  posv = first_rec[0]
  negv = first_rec[1]
  fid = first_rec[2]
  nrec = first_rec[3]
  ndims = first_rec[4]

  medv = (posv + negv) / 2.
  verp = posv - medv
  vern = negv - medv
  vers = (posv - negv) / 2.

  if ( (np.abs(medv-3000.) > 1.e-3) ):
    print('Something is wrong with the file %s' % fname)
    print('\tIs it a 3D file?')
    return np.nan

  # Stokes profiles:
  if (np.abs(fid-160904.) < 1.e-3):
    res = read_profile(fname, fmt_type=fmt_type)
  # Model atmosphere:
  elif (np.abs(fid-130904.) < 1.e-3):
    res = read_model(fname, fmt_type=fmt_type, devel=devel)
  elif (np.abs(fid-18060904.) < 1.e-3):
    res = read_rf(fname, itx=itx, fmt_type=fmt_type)
  else:
    print("\t")
    print("\t'firtez-dz.py':: read:")
    print("\tUnknown file type: %s" % (fname, ))
    return np.nan

  return res

def write_lsf(fname, data, fmt_type=np.float32):

  
  idd = 3330904
  v=3
  
  if (fmt_type != np.float32):
    print('Not implemented yet!')
    return np.nan
  
  # We write single precision numbers: 32 bits
  # This is 4 bytes

  nw=data.size
  nnum_to_write = 4. * nw
  nbyts_to_write = nnum_to_write * 4.
  
  # The maximum size of a fortran record is 2Gb
  nbyts_max_frec = 1. * 1024. * 1024. * 1024.
  nnum_max_frec = nbyts_max_frec / 4.
  
  # We set the maximum number of fortran rec to actual_max-10:
  nnum_rec = np.int(nnum_max_frec - 10.)
  nrecs_to_write = np.int32(np.ceil(nnum_to_write/nnum_rec))
  
  # FIRST, WE WRITE A FIRST RECORD WITH THE MANDATORY DATA:
  towrite = np.zeros(6, dtype=fmt_type)
  
  towrite[0] = 3000 + v
  towrite[1] = 3000 - v
  towrite[2] = idd * 1.
  towrite[3] = nrecs_to_write+1.
  towrite[4] = 1.  # NUMBER OF DIMENSIONS
  towrite[5] = nw * 1.
  
  #OPEN WRITING FILE:
  f=FortranFile(fname,'w')
  f.write_record(np.float32(towrite))
  #         SET OF DATA:
  offset = 0
  for it_nnn in range(nrecs_to_write):
    if (it_nnn != nrecs_to_write - 1):
      towrite = data[offset:offset+nnum_rec] * 1.
    else:
      towrite = data[offset:] * 1.
    f.write_record(np.float32(towrite))
    offset=offset+nnum_rec
  # CLOSE IT
  f.close()
  
  return


