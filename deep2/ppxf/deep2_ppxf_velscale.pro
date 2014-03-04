function deep2_ppxf_velscale
; jm09nov12ucsd - return the DEEP2 pixel size in km/s for performing
; the PPXF fitting; see J. Newman+13, Table 2

    light = 2.99792458D5        ; speed of light [km/s]
    fwhm2sig = 2D0*sqrt(2.0*alog(2.0))

    meanwave = 7800.0
    pixsize = 0.33D ; [mean pixel size, A]
    velscale = pixsize/meanwave*light

    velscale = 12.0             ; [km/s]

return, velscale
end
