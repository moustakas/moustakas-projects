function ediscs_ppxf_velscale, run34=run34
; jm10apr28ucsd - return the EDisCS pixel size in km/s for performing
; the PPXF fitting 
    light = 2.99792458D5        ; speed of light [km/s]
    meanwave = (8200.0-5200.0)/2.0+5200.0 ; mean wavelength
    if keyword_set(run34) then begin
       meanwave = mean([5120,8450])
       pixsize = 1.6D           ; [mean pixel size, A]
       velscale = pixsize/meanwave*light
       velscale = 70D           ; [km/s]
    endif else begin
       meanwave = mean([5300,8000])
       pixsize = 1.2D           ; [mean pixel size, A]
       velscale = pixsize/meanwave*light
       velscale = 55.0D         ; [km/s]
    endelse
return, velscale
end

