
FUNCTION dz_write_profiles, datastr, fname

IF (datastr.lb NE 160904) THEN BEGIN
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
plntot = LONG(datastr.nx)*LONG(datastr.ny)*LONG(datastr.nw)
lntot = 4l*plntot
to_be_written = lntot
loffs = lntot*0
nr = lntot/536870902
IF ((lntot MOD 536870902) NE 0) THEN nr = nr + 1
nr = nr + 3

data = REFORM([ $
      REFORM(TRANSPOSE(datastr.stki, [0,2,1]), 1,plntot) $
    , REFORM(TRANSPOSE(datastr.stkq, [0,2,1]), 1,plntot) $
    , REFORM(TRANSPOSE(datastr.stku, [0,2,1]), 1,plntot) $
    , REFORM(TRANSPOSE(datastr.stkv, [0,2,1]), 1,plntot)] $
    , lntot)


towrite = FLTARR(9)

towrite[0] = FLOAT(datastr.vp)
towrite[1] = FLOAT(datastr.vn)
towrite[2] = FLOAT(datastr.lb)
towrite[3] = FLOAT(nr)
towrite[4] = FLOAT(datastr.nd)
towrite[5] = FLOAT(datastr.nx)
towrite[6] = FLOAT(datastr.ny)
towrite[7] = FLOAT(datastr.nw)
towrite[8] = FLOAT(4)

OPENW, wunit, fname, /GET_LUN, /F77
WRITEU, wunit, towrite
WRITEU, wunit, datastr.index
WRITEU, wunit, datastr.waves

FOR i=0,nr-4 DO BEGIN
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
