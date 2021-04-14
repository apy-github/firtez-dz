
FUNCTION dz_write_models, datastr, fname

IF (datastr.lb NE 130904) THEN BEGIN
  PRINT, ''
  PRINT, ' Writing file: ', fname
  PRINT, ' Data label: ', datastr.lb
  PRINT, ' Is this a model atmosphere file?. Returning 0'
  PRINT, ''
  RETURN, 0
ENDIF

vv=(datastr.vp-datastr.vn)/2

IF (vv NE 3) THEN BEGIN

  PRINT, ''
  PRINT, ' Writing file: ', fname
  PRINT, ' Wrong model file version. Returning 0'
  PRINT, ''
  RETURN, 0

ENDIF

;  data:
plntot = LONG(datastr.nx)*LONG(datastr.ny)*LONG(datastr.nz)
lntot = LONG(datastr.npar)*plntot
to_be_written = lntot
loffs = lntot*0
nr = lntot/536870902
IF ((lntot MOD 536870902) NE 0) THEN nr = nr + 1
nr = nr + 1 ; Header

data = FLTARR(datastr.nz, datastr.ny, datastr.nx, datastr.npar)
data[*,*,*,0] = TRANSPOSE(REFORM(datastr.tem, datastr.nx, datastr.ny, datastr.nz), [2,1,0])
data[*,*,*,1] = TRANSPOSE(REFORM(datastr.pg, datastr.nx, datastr.ny, datastr.nz), [2,1,0])
data[*,*,*,2] = TRANSPOSE(REFORM(datastr.rho, datastr.nx, datastr.ny, datastr.nz), [2,1,0])
data[*,*,*,3] = TRANSPOSE(REFORM(datastr.bx, datastr.nx, datastr.ny, datastr.nz), [2,1,0])
data[*,*,*,4] = TRANSPOSE(REFORM(datastr.by, datastr.nx, datastr.ny, datastr.nz), [2,1,0])
data[*,*,*,5] = TRANSPOSE(REFORM(datastr.bz, datastr.nx, datastr.ny, datastr.nz), [2,1,0])
data[*,*,*,6] = TRANSPOSE(REFORM(datastr.vz, datastr.nx, datastr.ny, datastr.nz), [2,1,0])
data[*,*,*,7] = TRANSPOSE(REFORM(datastr.pel, datastr.nx, datastr.ny, datastr.nz), [2,1,0])
data[*,*,*,8] = TRANSPOSE(REFORM(datastr.mw, datastr.nx, datastr.ny, datastr.nz), [2,1,0])
data[*,*,*,9] = TRANSPOSE(REFORM(datastr.tau, datastr.nx, datastr.ny, datastr.nz), [2,1,0])
data[*,*,*,10] = TRANSPOSE(REFORM(datastr.x, datastr.nx, datastr.ny, datastr.nz), [2,1,0])
data[*,*,*,11] = TRANSPOSE(REFORM(datastr.y, datastr.nx, datastr.ny, datastr.nz), [2,1,0])
data[*,*,*,12] = TRANSPOSE(REFORM(datastr.z, datastr.nx, datastr.ny, datastr.nz), [2,1,0])

data = REFORM(data, lntot, /over)

towrite = FLTARR(9)

towrite[0] = FLOAT(datastr.vp)
towrite[1] = FLOAT(datastr.vn)
towrite[2] = FLOAT(datastr.lb)
towrite[3] = FLOAT(nr)
towrite[4] = FLOAT(datastr.nd)
towrite[5] = FLOAT(datastr.npar)
towrite[6] = FLOAT(datastr.nx)
towrite[7] = FLOAT(datastr.ny)
towrite[8] = FLOAT(datastr.nz)

OPENW, wunit, fname, /GET_LUN, /F77
WRITEU, wunit, towrite

FOR i=0,nr-2 DO BEGIN
  IF (to_be_written GT 536870902) THEN BEGIN
    lntot=536870902
  ENDIF ELSE BEGIN
    lntot=to_be_written
  ENDELSE
  WRITEU, wunit, data[loffs:loffs+lntot-1]
  loffs=loffs+lntot
  to_be_read=to_be_written-lntot
ENDFOR

FREE_LUN, wunit

RETURN,1

END
