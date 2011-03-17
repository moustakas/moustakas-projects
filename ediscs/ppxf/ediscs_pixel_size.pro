function ediscs_pixel_size, run34=run34
; jm10apr28ucsd - return the approximate linear pixel size for the
;   EDisCS spectra [A/pixel]
    if keyword_set(run34) then dwave = 1.6D else dwave = 1.2D
return, dwave
end
