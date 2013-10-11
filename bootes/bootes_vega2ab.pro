function bootes_vega2ab
; jm10jan27ucsd - return the BOOTES Vega-->AB conversion factors 
    filterlist = bootes_filterlist()
    vega2ab = k_vega2ab(filterlist=filterlist,/kurucz,/silent)

; LBC/u, LBC/Y, and zBootes magnitudes are already AB
    ok = where(strmatch(filterlist,'*lbc*',/fold) or $
      strmatch(filterlist,'*90prime_z*',/fold))
    vega2ab[ok] = 0.0

; use Reach et al. 2005 for the IRAC channels; see also
; http://ssc.spitzer.caltech.edu/irac/calib/
    irac = where(strmatch(filterlist,'*irac*',/fold))
    vega2ab[irac] = -2.5*alog10([280.9,179.7,115.0,64.13]*1D-23)-48.6

return, vega2ab
end
