function sg1120_filterlist, sdss=sdss
; jm09jul06nyu - build a QAplot for KCORRECT output; if PSFILE is
    if keyword_set(sdss) then begin
       filterlist = ['vimos_B','vimos_V','vimos_R',$
         'sdss_g0','sdss_r0','flamingos_Ks',$
         'sdss_u0','sdss_g0','sdss_r0','sdss_i0']+'.par'
    endif else begin
;      filterlist = ['vimos_B','vimos_V','vimos_R',$
;        'sdss_g0','sdss_r0']+'.par'
       filterlist = ['vimos_B','vimos_V','vimos_R',$
         'sdss_g0','sdss_r0','flamingos_Ks']+'.par'
    endelse
return, filterlist
end
