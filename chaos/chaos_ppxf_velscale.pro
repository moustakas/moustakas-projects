function chaos_ppxf_velscale
; jm13apr01siena - return the CHAOS pixel size in km/s for performing
; the PPXF fitting 
    light = 2.99792458D5        ; speed of light [km/s]
    fwhm2sig = 2D0*sqrt(2.0*alog(2.0))
    meanwave = (9200.0-3700.0)/2.0+3700.0 ; CHAOS
;   fwhm = 6.0 ; [FWHM instrumental resolution, A]
;   velscale = fwhm/meanwave/fwhm2sig*light
;   velscale = 120.0 ; [km/s]
    pixsize = 1.2D ; [mean pixel size, A]
    velscale = pixsize/meanwave*light
    velscale = 55.0D ; [km/s]
return, velscale
end
