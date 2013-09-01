function redbaryons_zbins, nzbins, zmin=zmin, zmax=zmax, $
  area=area
; jm13aug27siena - preselected redshift bins
; redshift intervals correspond to ~0.62 Gyr of cosmic time 

;; figure out the optimal redshift bins; this
;    zmin = 0.1
;    zmax = 0.5
;    nz = 6.0
;    dage = (getage(zmin)-getage(zmax))/nz
;    niceprint, reverse(getredshift(findgen(nz)*dage+dage+getage(zmax))), $
;      reverse(getredshift(findgen(nz)*dage+getage(zmax)))
;      0.10000000        0.15340778
;      0.15340778        0.21107179
;      0.21107179        0.27364513
;      0.27364513        0.34195873
;      0.34195873        0.41695154
;      0.41695154        0.50001603

    zzmin = [0.10,0.15,0.20,0.25,0.30,0.35]
    zzmax = [0.15,0.20,0.25,0.30,0.35,0.40]
;   splog, mean(getage(zzmin)-getage(zzmax))
    
    zmin = min(zzmin)
    zmax = max(zzmax)
    nzbins = n_elements(zzmin)
    zbins = replicate({zbin: 0.0, zlo: 0.0, zup: 0.0, $
      zsigma: 0.0, time: 0.0, vol: 0.0},nzbins)
    zbins.zlo = zzmin
    zbins.zup = zzmax
    zbins.zbin = (zzmax-zzmin)/2.0+zzmin
    zbins.zsigma = (zbins.zup-zbins.zlo)/2.0
    zbins.time = getage(zzmin)-getage(zzmax)

    h100 = 0.7
    area = 10400.0*(!pi/180.0)^2
    for iz = 0, nzbins-1 do zbins[iz].vol = (area/3.0)*$
      (lf_comvol(zbins[iz].zup)-lf_comvol(zbins[iz].zlo))*(1.0/h100)^3.0 ; h=1-->h=0.7
    
return, zbins
end
    
