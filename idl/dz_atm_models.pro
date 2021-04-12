
FUNCTION dz_atm_models, tem=tem, pg=pg, rho=rho,  bx=bx,  by=by,  bz=bz $
    ,  vz=vz,  tau=tau,  pel=pel,  mw=mw,  x=x,  y=y,  z=z, nx=nx, ny=ny, nz=nz
 
dims = 0
tofill = intarr(13) + 1
IF (N_ELEMENTS(tem) NE 0) THEN BEGIN
  dims = SIZE(tem, /dim)
  tofill[0] = 0
ENDIF
IF (N_ELEMENTS(pg) NE 0) THEN BEGIN
  IF (dims EQ 0) THEN dims = SIZE(pg, /dim)
  tofill[1] = 0
ENDIF
IF (N_ELEMENTS(rho) NE 0) THEN BEGIN
  IF (dims EQ 0) THEN dims = SIZE(rho, /dim)
  tofill[2] = 0
ENDIF
IF (N_ELEMENTS(bx) NE 0) THEN BEGIN
  IF (dims EQ 0) THEN dims = SIZE(bx, /dim)
  tofill[3] = 0
ENDIF
IF (N_ELEMENTS(by) NE 0) THEN BEGIN
  IF (dims EQ 0) THEN dims = SIZE(by, /dim)
  tofill[4] = 0
ENDIF
IF (N_ELEMENTS(bz) NE 0) THEN BEGIN
  IF (dims EQ 0) THEN dims = SIZE(bz, /dim)
  tofill[5] = 0
ENDIF
IF (N_ELEMENTS(vz) NE 0) THEN BEGIN
  IF (dims EQ 0) THEN dims = SIZE(vz, /dim)
  tofill[6] = 0
ENDIF
IF (N_ELEMENTS(tau) NE 0) THEN BEGIN
  IF (dims EQ 0) THEN dims = SIZE(tau, /dim)
  tofill[7] = 0
ENDIF
IF (N_ELEMENTS(pel) NE 0) THEN BEGIN
  IF (dims EQ 0) THEN dims = SIZE(pel, /dim)
  tofill[8] = 0
ENDIF
IF (N_ELEMENTS(mw) NE 0) THEN BEGIN
  IF (dims EQ 0) THEN dims = SIZE(mw, /dim)
  tofill[9] = 0
ENDIF
IF (N_ELEMENTS(x) NE 0) THEN BEGIN
  IF (dims EQ 0) THEN dims = SIZE(x, /dim)
  tofill[10] = 0
ENDIF
IF (N_ELEMENTS(y) NE 0) THEN BEGIN
  IF (dims EQ 0) THEN dims = SIZE(y, /dim)
  tofill[11] = 0
ENDIF
IF (N_ELEMENTS(z) NE 0) THEN BEGIN
  IF (dims EQ 0) THEN dims = SIZE(z, /dim)
  tofill[12] = 0
ENDIF

IF (TOTAL(1-tofill) EQ 0) THEN BEGIN
  IF ( (N_ELEMENTS(nx) EQ 0) OR (N_ELEMENTS(ny) EQ 0) OR (N_ELEMENTS(nz) EQ 0) ) THEN BEGIN
    PRINT, ''
    PRINT, ''
    PRINT, ''
    RETURN, 1./0.
  ENDIF ELSE BEGIN
    dims = [nx, ny, nz]
  ENDELSE
ENDIF

FOR i=0,N_ELEMENTS(tofill)-1 DO BEGIN
  IF (tofill[i] EQ 1) THEN BEGIN
    CASE i of
      0: tem = fltarr(dims)
      1: pg = fltarr(dims)
      2: rho = fltarr(dims)
      3: bx = fltarr(dims)
      4: by = fltarr(dims)
      5: bz = fltarr(dims)
      6: vz = fltarr(dims)
      7: tau = fltarr(dims)
      8: pel = fltarr(dims)
      9: mw = fltarr(dims)
      10: x = fltarr(dims)
      11: y = fltarr(dims)
      12: z = fltarr(dims)
    ENDCASE
  ENDIF
ENDFOR

nx = dims[0]
ny = dims[1]
nz = dims[2]

RETURN, {name:'' $
    , vp:3003 $
    , vn:2997 $
    , lb:130904 $
    , nd:LONG(4) $
    , npar:LONG(13) $
    , nx:LONG(nx) $
    , ny:LONG(ny) $
    , nz:LONG(nz) $
    , tem: REFORM(FLOAT(tem), nx, ny, nz)$
    , pg: REFORM(FLOAT(pg), nx, ny, nz)$
    , rho: REFORM(FLOAT(rho), nx, ny, nz)$
    , bx: REFORM(FLOAT(bx), nx, ny, nz)$
    , by: REFORM(FLOAT(by), nx, ny, nz)$
    , bz: REFORM(FLOAT(bz), nx, ny, nz)$
    , vz: REFORM(FLOAT(vz), nx, ny, nz)$
    , tau: REFORM(FLOAT(tau), nx, ny, nz)$
    , pel: REFORM(FLOAT(pel), nx, ny, nz)$
    , mw: REFORM(FLOAT(mw), nx, ny, nz)$
    , x: REFORM(FLOAT(x), nx, ny, nz)$
    , y: REFORM(FLOAT(y), nx, ny, nz)$
    , z: REFORM(FLOAT(z), nx, ny, nz)$
    }

END
