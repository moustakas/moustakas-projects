function wise_ppxf_velscale
; jm11apr13ucsd - return the SDSS pixel size in km/s for performing
; the PPXF fitting 
    velscale = 1D-4*im_light(/km)*alog(10)
return, velscale
end
