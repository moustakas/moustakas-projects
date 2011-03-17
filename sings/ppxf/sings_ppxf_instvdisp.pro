function sings_ppxf_instvdisp, bc03_instvdisp=bc03_instvdisp
; jm09nov12ucsd - convolve the BC03 model templates by the SINGS
; instrumental velocity dispersion (see BUILD_SINGS_PPXF_TEMPLATES);
; assume to be constant with wavelength
    light = 2.99792458D5        ; speed of light [km/s]
    fwhm2sig = 2D0*sqrt(2.0*alog(2.0))

    meanwave = (7000.0-3600.0)/2.0+3600.0 ; mean wavelength
    fwhm = 8.5 ; [FWHM instrumental resolution, A]
    instvdisp = fwhm/meanwave/fwhm2sig*light
    instvdisp = 200.0 ; [km/s]
;   instvdisp = 100.0 ; [km/s]

    if arg_present(bc03_instvdisp) then begin
       bc03_fwhm = 3.0 ; [A]
       bc03_instvdisp = 3.0/meanwave/fwhm2sig*light
       bc03_instvdisp = 70.0 ; [km/s]
    endif

return, instvdisp
end
