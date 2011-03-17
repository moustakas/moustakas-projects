function cosmicimf_zbins, nzbins, lf24=lf24
; jm10mar17ucsd - redshift bins for the MFs, by default, or for the
; L(24) LFs if /lf24

; choose redshifts bins corresponding to equal intervals of cosmic
; time for both the stellar mass and SFR functions:
;   niceprint, dtime(0.05,0.75,6)
    if keyword_set(lf24) then begin
       nzbins = 4
       zbins = replicate({zbin: 0.0, zbin_err: 0.0, zlo: 0.0, zup: 0.0},nzbins)
       zbins.zlo = [0.05,0.20,0.35,0.50]
       zbins.zup = [0.20,0.35,0.50,0.65]
;      zbins.zlo = [0.05,0.15,0.25,0.35,0.45,0.55]
;      zbins.zup = [0.15,0.25,0.35,0.45,0.55,0.75]
;      zbins.zlo = [0.050,0.130,0.220,0.323,0.443,0.584]
;      zbins.zup = [0.130,0.220,0.323,0.443,0.584,0.753]
;      zbins.zlo = [0.147,0.392,0.220]
;      zbins.zup = [0.260,0.554,0.443]
    endif else begin
       nzbins = 6
       zbins = replicate({zbin: 0.0, zbin_err: 0.0, zlo: 0.0, zup: 0.0},nzbins)
       zbins.zlo = [0.05,0.15,0.25,0.35,0.45,0.55]
       zbins.zup = [0.15,0.25,0.35,0.45,0.55,0.75]
;      zbins.zlo = [0.050,0.130,0.220,0.323,0.443,0.584]
;      zbins.zup = [0.130,0.220,0.323,0.443,0.584,0.753]
;      zbins.zlo = [0.05,0.260,0.554,0.05,0.443]
;      zbins.zup = [0.147,0.392,0.753,0.220,0.753]
    endelse
    zbins.zbin = total([[zbins.zlo],[zbins.zup]],2)/2.0
    zbins.zbin_err = (zbins.zup-zbins.zlo)/2.0
return, zbins
end
