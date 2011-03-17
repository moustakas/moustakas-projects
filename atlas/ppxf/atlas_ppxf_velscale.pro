function atlas_ppxf_velscale
; jm09dec17ucsd - return the ATLAS pixel size in km/s for performing
; the PPXF fitting 
    light = 2.99792458D5        ; speed of light [km/s]
    meanwave = (7000.0-3600.0)/2.0+3600.0 ; mean wavelength
    pixsize = 2.75D ; [mean pixel size, A]
    velscale = pixsize/meanwave*light
    velscale = 155.0D ; [km/s]
return, velscale
end
