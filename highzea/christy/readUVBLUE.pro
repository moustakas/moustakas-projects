function readUVBLUE, filename0, VERBOSE=verbose

;; This IDL function reads the file of a high-resolution spectrum of the
;; UBLUE library and stores the data in the structure s, which
;; contains the effective temperature (Teff), the surface gravity log(g)
;; (grav), and the metallicity [M/H] (metal) of the model, and the spectral 
;; resolution (resolu), the wavelength (lam), the normalized flux (fnlam), the
;; residual flux (res), the emergent flux (flam), and the continuum flux
;; (fclam) of the spectrum.
;; Wavelengths are in angstrom.
;; The emergent flux is in erg s^-1 cm^-2 A^-1
;; The residual flux is set to zero where the continuum flux is zero.
;; The bolometric integral of the normalized flux is 1.
;; The files to be read can be compressed (must have extension .gz or .bz2)
;;
;; written by Emanuele Bertone, INAOE
;; version: April 1, 2005


if N_PARAMS() ne 1 then begin
    PRINT,'Syntax - s = readUVBLUE(filename [,/VERBOSE])'
    PRINT,'   e.g. s = readUVBLUE(''./t06500g50p00k2.flx.gz'')'
    return, -1
endif

lungh = STRLEN(filename0)
if STRMID(filename0,lungh-2,2) eq 'gz' then begin
    if KEYWORD_SET( VERBOSE ) then PRINT, 'COMPRESSED FILE.  Uncompressing...'
    spawn, STRING('gunzip -c '+filename0+' > ./tmp_spectrum'),ierr
    filename = './tmp_spectrum'
endif else if STRMID(filename0,lungh-3,3) eq 'bz2' then begin
    if KEYWORD_SET( VERBOSE ) then PRINT, 'COMPRESSED FILE.  Uncompressing...'
    spawn, STRING('bunzip2 -c '+filename0+' > ./tmp_spectrum'),ierr
    filename = './tmp_spectrum'
endif else filename = filename0

OPENR, unit, filename, /get_lun

teff   = 0.
grav   = 0.
metal  = 0.
wlbwg  = 0.d0
wlend  = 0.d0
resolu = 0.
ratio  = 0.d0
nwl    = 0L
w0 = 0.d0
f0 = 0.d0
r0 = 0.d0

READF, unit, teff, grav, metal
READF, unit, nwl, wlbeg, wlend
READF, unit, resolu, ratio

if KEYWORD_SET( VERBOSE ) then PRINT, 'Reading spectrum... Teff:', teff, '  Log g:', grav, '  [M/H]:', metal, format='(a, f8.0, a, f4.1, a, f5.1)'

s = {teff: 0., grav: 0., metal: 0., resolu: 0., lam: DBLARR(nwl), fnlam: DBLARR(nwl), res: DBLARR(nwl), flam: DBLARR(nwl), fclam: DBLARR(nwl)}

s.teff   = teff
s.grav   = grav
s.metal  = metal
s.resolu = resolu
s.lam = DOUBLE(wlbeg) * DOUBLE(ratio)^DINDGEN(nwl)

data = FLTARR(2,nwl)
READF, unit, data
s.flam  = data[0,*]
s.fclam = data[1,*]

sigma = 5.6697d-5             ; Stefan-Boltzmann constant
                              ; value used by Kurucz in his codes

s.fnlam = s.flam / (sigma * s.teff^4)
ind0 = WHERE(s.fclam eq 0., complement=ind1, cc)
if cc ge 1 then s.res[ind0] = 0.
s.res[ind1]   = s.flam[ind1] / s.fclam[ind1]

CLOSE, unit
FREE_LUN, unit

if STRMID(filename0,lungh-2,2) eq 'gz' or STRMID(filename0,lungh-3,3) eq 'bz2' then SPAWN,'/bin/rm -f ./tmp_spectrum',ierr

if KEYWORD_SET( VERBOSE ) then PRINT, nwl, ' wavelength points read.'

return, s

end
