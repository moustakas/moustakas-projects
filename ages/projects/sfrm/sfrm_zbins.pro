function sfrm_zbins, nzbins
; jm10feb05ucsd
    nzbins = 6
    zbins = replicate({zbin: 0.0, zlo: 0.0, zup: 0.0},nzbins)
    zbins.zlo = [0.05,0.15,0.25,0.35,0.45,0.55]
    zbins.zup = [0.15,0.25,0.35,0.45,0.55,0.75]
    zbins.zbin = [0.1,0.2,0.3,0.4,0.5,0.65]
return, zbins
end
    
