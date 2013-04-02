function chaos_ppxf_instvdisp, bc03_instvdisp=bc03_instvdisp
; jm13apr01siena - convolve the BC03 model templates by the CHAOS
; instrumental velocity dispersion (see BUILD_CHAOS_PPXF_TEMPLATES);
; assume to be constant with wavelength
    light = 2.99792458D5        ; speed of light [km/s]
    fwhm2sig = 2D0*sqrt(2.0*alog(2.0))
    meanwave = (9200.0-3700.0)/2.0+3700.0 ; AGES
    fwhm = 6.0 ; [FWHM instrumental resolution, A]
    instvdisp = fwhm/meanwave/fwhm2sig*light
    instvdisp = 80.0 ; [km/s]
;   instvdisp = 120.0 ; [km/s]

    if arg_present(bc03_instvdisp) then begin
       bc03_fwhm = 3.0 ; [A]
       bc03_instvdisp = 3.0/meanwave/fwhm2sig*light
       bc03_instvdisp = 60.0 ; [km/s]
    endif

return, instvdisp
end
