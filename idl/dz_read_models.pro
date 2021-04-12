
FUNCTION dz_read_models, fname

hdr1 = FLTARR(5)
OPENR, runit, fname, /GET_LUN, /F77
READU, runit, hdr1
FREE_LUN, runit

vp=LONG(hdr1[0])
vn=LONG(hdr1[1])
lb=LONG(hdr1[2])
nr=LONG(hdr1[3])
nd=LONG(hdr1[4])

IF (lb NE 130904) THEN BEGIN
  PRINT, ''
  PRINT, ' Reading file: ', fname
  PRINT, ' Is this a model atmosphere file?. Returning 0'
  PRINT, ''
  RETURN, 0
ENDIF

vv=(vp-vn)/2

IF (vv NE 3) THEN BEGIN

  PRINT, ''
  PRINT, ' Reading file: ', fname
  PRINT, ' Wrong model file version. Returning 0'
  PRINT, ''
  RETURN, 0

ENDIF


hdr1 = FLTARR(5+nd)
OPENR, runit, fname, /GET_LUN, /F77
READU, runit, hdr1

npar=LONG(hdr1[5])
nx=LONG(hdr1[6])
ny=LONG(hdr1[7])
nz=LONG(hdr1[8])

lntot=LONG(npar)*LONG(nx)*LONG(ny)*LONG(nz)
to_be_read=lntot
loffs=lntot*0
data = FLTARR(lntot)

FOR i=1,nr-1 DO BEGIN
  IF (to_be_read GT 536870902) THEN BEGIN
    lntot=536870902
  ENDIF ELSE BEGIN
    lntot=to_be_read
  ENDELSE
  it_data = FLTARR(lntot)
  READU, runit, it_data
  data[loffs:loffs+lntot-1]=it_data
  loffs=loffs+lntot
  to_be_read=to_be_read-lntot
ENDFOR
FREE_LUN, runit

data = REFORM(data, nz,ny,nx,npar, /over)

RETURN, {name:'' $
    , vn:vn $
    , vp:vp $
    , lb:lb $
    , nd:nd $
    , nx:nx $
    , ny:ny $
    , nz:nz $
    , npar:npar $
    , tem:REFORM(TRANSPOSE(REFORM(data[*,*,*,0],nz,ny,nx),[2,1,0]),nx,ny,nz) $
    , pg:REFORM(TRANSPOSE(REFORM(data[*,*,*,1],nz,ny,nx),[2,1,0]),nx,ny,nz) $
    , rho:REFORM(TRANSPOSE(REFORM(data[*,*,*,2],nz,ny,nx),[2,1,0]),nx,ny,nz) $
    , bx:REFORM(TRANSPOSE(REFORM(data[*,*,*,3],nz,ny,nx),[2,1,0]),nx,ny,nz) $
    , by:REFORM(TRANSPOSE(REFORM(data[*,*,*,4],nz,ny,nx),[2,1,0]),nx,ny,nz) $
    , bz:REFORM(TRANSPOSE(REFORM(data[*,*,*,5],nz,ny,nx),[2,1,0]),nx,ny,nz) $
    , vz:REFORM(TRANSPOSE(REFORM(data[*,*,*,6],nz,ny,nx),[2,1,0]),nx,ny,nz) $
    , tau:REFORM(TRANSPOSE(REFORM(data[*,*,*,7],nz,ny,nx),[2,1,0]),nx,ny,nz) $
    , pel:REFORM(TRANSPOSE(REFORM(data[*,*,*,8],nz,ny,nx),[2,1,0]),nx,ny,nz) $
    , mw:REFORM(TRANSPOSE(REFORM(data[*,*,*,9],nz,ny,nx),[2,1,0]),nx,ny,nz) $
    , x:REFORM(TRANSPOSE(REFORM(data[*,*,*,10],nz,ny,nx),[2,1,0]),nx,ny,nz) $
    , y:REFORM(TRANSPOSE(REFORM(data[*,*,*,11],nz,ny,nx),[2,1,0]),nx,ny,nz) $
    , z:REFORM(TRANSPOSE(REFORM(data[*,*,*,12],nz,ny,nx),[2,1,0]),nx,ny,nz) $
    }
END
