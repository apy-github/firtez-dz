
FUNCTION dz_stk_profiles, stki, stkq, stku, stkv

dims = SIZE(stki,/dim)

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
    , stki: FLOAT(stki)$
    , stkq: FLOAT(stkq)$
    , stku: FLOAT(stku)$
    , stkv: FLOAT(stkv)$
    }

END
