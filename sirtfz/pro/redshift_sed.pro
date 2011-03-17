pro redshift_sed, z, wave_0, lum_0, wave_z, flux_z, jansky=jansky
;+
; NAME:
;	REDSHIFT_SED
;
; PURPOSE:
;	Redshift a rest-frame SED (in W/Hz) to any redshift.
;
; CALLING SEQUENCE:
;	redshift_sed, z, wave_0, lum_0, wave_z, flux_z, mjansky=mjansky
;
; INPUTS:
;	z       - redshift
;	wave_0  - rest-frame wavelength (any units)
;	lum_0   - monochromatic luminosity (W/Hz)
;
; KEYWORD PARAMETERS:
;	jansky - return the flux density in Jy
;
; OUTPUTS:
;	wave_z  - redshifted wavelength vector
;	flux_z  - flux density in W/m^2/Hz (or Jy)
;
; COMMON BLOCKS:
;	The common block COSMOLOGY is implicit via the DLUMINOSITY()
;	routine.  
;
; COMMENTS:
;	See Madau 1995, ApJ, 441, 18, Eq. (3) for the relation between
;	mean specific flux and redshift.
;
; PROCEDURES USED:
;	DLUMINOSITY()
;
; MODIFICATION HISTORY:
;	J. Moustakas, 2001 Sep 26, U of A, written
;	jm01jul30uofa, documented
;-

    if n_params() eq 0L then begin
       print, 'Syntax - redshift_sed, z, wave_0, lum_0, wave_z, flux_z, mjansky=mjansky'
       return
    endif

    if (z le 0.0) then begin
       print, 'Setting the redshift to zero.'
       z = 0.0
    endif

    wave_z = wave_0 * (1.0+z)
    dlum = dluminosity(z) * 3.085678D16                  ; luminosity distance (m)
    flux_z = (1.0+z) * lum_0 / (4.0D*!dpi*dlum*dlum)     ; flux density (W/m^2/Hz)
    if keyword_set(jansky) then flux_z = flux_z / 1D-26 ; flux in Jy

return
end
