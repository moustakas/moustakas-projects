function mz_zbins, nzbins, zmin=zmin, zmax=zmax, sdss=sdss
; jm10may08ucsd
    if keyword_set(sdss) then begin
       nzbins = 1
       zbins = {zbin: 0.0, zlo: 0.0, zup: 0.0}
       zbins.zlo = 0.033
       zbins.zup = 0.25
       zbins.zbin = 0.1415
    endif else begin
       nzbins = 6
       zbins = replicate({zbin: 0.0, zlo: 0.0, zup: 0.0},nzbins)
       zbins.zlo = [0.05,0.15,0.25,0.35,0.45,0.55]
       zbins.zup = [0.15,0.25,0.35,0.45,0.55,0.75]
       zbins.zbin = [0.1,0.2,0.3,0.4,0.5,0.65]
    endelse
    zmin = min(zbins.zlo)
    zmax = max(zbins.zup)
       
return, zbins
end
    
