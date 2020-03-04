
FUNCTION dz_read_profiles, fname

hdr1 = FLTARR(5)
OPENR, runit, fname, /GET_LUN, /F77
READU, runit, hdr1
FREE_LUN, runit

vp=LONG(hdr1[0])
vn=LONG(hdr1[1])
lb=LONG(hdr1[2])
nr=LONG(hdr1[3])
nd=LONG(hdr1[4])

IF (lb NE 160904) THEN BEGIN
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

nx=LONG(hdr1[5])
ny=LONG(hdr1[6])
nw=LONG(hdr1[7])
ns=LONG(hdr1[8])

wave = FLTARR(nw)
indx = FLTARR(nw)
lntot=LONG(ns)*LONG(nx)*LONG(ny)*LONG(nw)
to_be_read=lntot
loffs=lntot*0
data = FLTARR(lntot)

FOR i=1,nr-1 DO BEGIN
  IF (i EQ 1) THEN BEGIN
    READU, runit, indx
  ENDIF ELSE IF (i EQ 2) THEN BEGIN
    READU, runit, wave
  ENDIF ELSE BEGIN
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
  ENDELSE
ENDFOR
FREE_LUN, runit

data = REFORM(data, ns,nw,ny,nx, /over)

RETURN, {name:'' $
    , vn:vn $
    , vp:vp $
    , lb:lb $
    , nd:nd $
    , nw:nw $
    , nx:nx $
    , ny:ny $
    , waves:wave $
    , index:indx $
    , stki:REFORM(TRANSPOSE(REFORM(data[0,*,*,*],nw,ny,nx),[0,2,1]),nw,nx,ny) $
    , stkq:REFORM(TRANSPOSE(REFORM(data[1,*,*,*],nw,ny,nx),[0,2,1]),nw,nx,ny) $
    , stku:REFORM(TRANSPOSE(REFORM(data[2,*,*,*],nw,ny,nx),[0,2,1]),nw,nx,ny) $
    , stkv:REFORM(TRANSPOSE(REFORM(data[3,*,*,*],nw,ny,nx),[0,2,1]),nw,nx,ny) $
    }

END
