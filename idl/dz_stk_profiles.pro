
FUNCTION dz_stk_profiles, stki=stki, stkq=stkq, stku=stku, stkv=stkv, nx=nx, ny=ny, nw=nw
 
dims = 0
tofill = intarr(4) + 1
IF (N_ELEMENTS(stki) NE 0) THEN BEGIN
  dims = SIZE(stki, /dim)
  tofill[0] = 0
ENDIF
IF (N_ELEMENTS(stkq) NE 0) THEN BEGIN
  IF (dims EQ 0) THEN dims = SIZE(stkq, /dim)
  tofill[1] = 0
ENDIF
IF (N_ELEMENTS(stku) NE 0) THEN BEGIN
  IF (dims EQ 0) THEN dims = SIZE(stku, /dim)
  tofill[2] = 0
ENDIF
IF (N_ELEMENTS(stkv) NE 0) THEN BEGIN
  IF (dims EQ 0) THEN dims = SIZE(stkv, /dim)
  tofill[3] = 0
ENDIF

IF (TOTAL(1-tofill) EQ 0) THEN BEGIN
  IF ( (N_ELEMENTS(nx) EQ 0) OR (N_ELEMENTS(ny) EQ 0) OR (N_ELEMENTS(nw) EQ 0) ) THEN BEGIN
    PRINT, ''
    PRINT, ''
    PRINT, ''
    RETURN, 1./0.
  ENDIF ELSE BEGIN
    dims = [nw, nx, ny]
  ENDELSE
ENDIF

FOR i=0,N_ELEMENTS(tofill)-1 DO BEGIN
  IF (tofill[i] EQ 1) THEN BEGIN
    CASE i of
      0: stki = fltarr(dims)
      1: stkq = fltarr(dims)
      2: stku = fltarr(dims)
      3: stkv = fltarr(dims)
    ENDCASE
  ENDIF
ENDFOR

nw = dims[0]
nx = dims[1]
ny = dims[2]

RETURN, {name:'' $
    , vp:3003 $
    , vn:2997 $
    , lb:160904 $
    , nd:LONG(4) $
    , nw:LONG(nw) $
    , nx:LONG(nx) $
    , ny:LONG(ny) $
    , waves:FLTARR(nw) $
    , index:FLTARR(nw) $
    , stki: REFORM(FLOAT(stki), nw, nx, ny)$
    , stkq: REFORM(FLOAT(stkq), nw, nx, ny)$
    , stku: REFORM(FLOAT(stku), nw, nx, ny)$
    , stkv: REFORM(FLOAT(stkv), nw, nx, ny)$
    }

END
