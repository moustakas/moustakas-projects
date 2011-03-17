function ediscs_ppxf_instvdisp, bc03_instvdisp=bc03_instvdisp, run34=run34
; jm10apr28ucsd - convolve the BC03 model templates by the EDISCS
; instrumental velocity dispersion (see BUILD_EDISCS_PPXF_TEMPLATES);
; assume to be constant with wavelength
    light = 2.99792458D5        ; speed of light [km/s]
    fwhm2sig = 2D0*sqrt(2.0*alog(2.0))

    fwhm = 6.0 ; [FWHM instrumental resolution, A]
    if keyword_set(run34) then begin
       meanwave = mean([5120,8450])
       instvdisp = fwhm/meanwave/fwhm2sig*light
       instvdisp = 115.0 ; [km/s]
    endif else begin
       meanwave = mean([5300,8000])
       instvdisp = fwhm/meanwave/fwhm2sig*light
       instvdisp = 115.0 ; [km/s]
    endelse    
    
    if arg_present(bc03_instvdisp) then begin
       bc03_fwhm = 3.0 ; [A]
       bc03_instvdisp = 3.0/meanwave/fwhm2sig*light
       bc03_instvdisp = 70.0 ; [km/s]
    endif

return, instvdisp
end
