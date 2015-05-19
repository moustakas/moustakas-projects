function desi_templates_ppxf_velscale
; jm09dec17ucsd - return the ATLAS pixel size in km/s for performing
; the PPXF fitting 
    light = 2.99792458D5        ; speed of light [km/s]
    meanwave = (9500.0-3650.0)/2.0+3650.0 ; mean wavelength
    pixsize = 0.5D ; [mean pixel size, A]
    velscale = pixsize/meanwave*light
    velscale = 25D ; [km/s]
return, velscale
end
