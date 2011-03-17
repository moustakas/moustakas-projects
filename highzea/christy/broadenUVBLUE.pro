function broadenUVBLUE, lam0, flam0, fwhm=fwhm, resolu=resolu, vel=vel, nsigma=nsigma, inputresolu=inputresolu, inputreslam=inputreslam

;; This IDL function broadens a UVBLUE spectrum to the desired spectral
;; resolution by convolution with a Gaussian kernel. The output spectral
;; resolution can be expressed as a FWHM, an inverse resolution R, or a
;; velocity value.
;; The input data are the array containing the wavelength points and the array
;; containing the flux value. The arrays must have the same dimension.
;;
;; written by Emanuele Bertone, INAOE
;; version: March 25, 2005

ON_ERROR,2                      ;Return to caller

nkeyw = 0
if KEYWORD_SET( FWHM )   then nkeyw = nkeyw + 1
if KEYWORD_SET( RESOLU ) then nkeyw = nkeyw + 1
if KEYWORD_SET( VEL )    then nkeyw = nkeyw + 1

if N_PARAMS() ne 2 and nkeyw lt 1 then begin
    PRINT,'Syntax - fbroad = broadenUVBLUE(wavelength array, flux array, FWHM=fwhm, RESOLU=resolu, VEL=vel [, INPUTRESOLU=inputresolu, INPUTRESLAM=inputreslam, NSIGMA=nsigma])'
    PRINT, '  Wavelengths must be in angstrom'
    PRINT, '  RESOLU: output resolution [R = lambda/(Delta lambda)]'
    PRINT, '  FWHM:   output fwhm [in angstrom] (exact only in the interval center)'
    PRINT, '  VEL:    output velocity dispersion [in km/s]'
    PRINT, '  INPUTRESOLU: spectral resolution of the input spectrum'
    PRINT, '  INPUTRESLAM: resolution corresponding to the wavelength step'
    PRINT, '  NSIGMA: width of the gaussian kernel in sigma units; default is 5.'
    PRINT, ' Notes:'
    PRINT, '  Only one of FWHM, RESOLU, VEL keyword must be present'
    PRINT, '  Convolution kernel is Gaussian'
    PRINT, '  Wavelength step must be proportional to lam/inputreslam'
    PRINT, '  Default values are: inputresolu=10000., inputreslam=40000., nsigma=5'
    return, -1
endif

if not KEYWORD_SET( INPUTRESOLU ) then inputresolu = 10000.d0
if not KEYWORD_SET( INPUTRESLAM ) then inputreslam = 40000.d0
if not KEYWORD_SET( NSIGMA ) then nsigma = 5.

c = 2.99792458d5

lam  = DOUBLE(lam0[*])
flam = DOUBLE(flam0[*])

if KEYWORD_SET( FWHM ) then begin
;; if this keyword is set the output FWHM of the broadened spectrum is exact
;; only in the center of the chose wavelength interval
    sigma = DOUBLE(fwhm/2.3548)
    nwl = N_ELEMENTS(lam)
    dlam = lam[LONG(nwl/2L)] - lam[LONG(nwl/2L)-1L]
    nk = LONG(2.0*nsigma*sigma/dlam)
;; check and set to odd the wavelength points of the kernel
    resto = nk mod 2
    imed = nk/2
    if resto eq 0 then nk = nk + 1
    kern = FLTARR(nk)
    fimed = DOUBLE(imed)
    lamk = (FINDGEN(nk) - fimed)*dlam
    kern = EXP(-0.5*(lamk/sigma)^2)
    somma = TOTAL(kern)
    kern = kern/somma
endif else if KEYWORD_SET( RESOLU ) then begin
    resolu2 = SQRT((DOUBLE(resolu)*inputresolu)^2/(DOUBLE(inputresolu)^2 - resolu^2))
    imax = FIX(nsigma/2.3548*DOUBLE(inputreslam/resolu2))
    ik = DINDGEN(2*imax + 1) - imax
    kern = EXP(-0.5*(2.3548*ik*DOUBLE(resolu2/inputreslam))^2)
    somma = TOTAL(kern)
    kern = kern/somma
endif else if KEYWORD_SET( VEL ) then begin
    resolu = c / DOUBLE(vel)
    resolu2 = SQRT((DOUBLE(resolu)*inputresolu)^2/(DOUBLE(inputresolu)^2 - resolu^2))
    imax = FIX(nsigma/2.3548*DOUBLE(inputreslam/resolu2))
    ik = DINDGEN(2*imax + 1) - imax
    kern = EXP(-0.5*(2.3548*ik*DOUBLE(resolu2/inputreslam))^2)
    somma = TOTAL(kern)
    kern = kern/somma
endif
broad = CONVOL(flam, kern, /edge_truncate)

return, broad

end
