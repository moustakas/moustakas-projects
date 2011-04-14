function wise_ppxf_instvdisp, bc03_instvdisp=bc03_instvdisp
; jm11apr13ucsd - convolve the BC03 model templates by the SDSS
; instrumental velocity dispersion (see BUILD_WISE_PPXF_TEMPLATES);
; assume to be constant with wavelength; see Tremonti+04 for details 
    light = 2.99792458D5        ; speed of light [km/s]
    fwhm2sig = 2D0*sqrt(2.0*alog(2.0))

    meanwave = 5000.0
    fwhm = 2.4 ; [FWHM instrumental resolution, A] @5000 a
    instvdisp = fwhm/meanwave/fwhm2sig*light
    instvdisp = 65.0 ; [km/s]

    if arg_present(bc03_instvdisp) then begin
       bc03_fwhm = 3.0 ; [A]
       bc03_instvdisp = 3.0/meanwave/fwhm2sig*light
       bc03_instvdisp = 60.0 ; [km/s]
    endif

return, instvdisp
end
