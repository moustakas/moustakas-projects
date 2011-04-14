function mzbpt_upper_limits, iratios1, mz_ispec, doplot=doplot
; internal support routine for BUILD_MZ_EMLINE; classify using the BPT
; diagram and upper limits

; deal with upper limits; consider three combinations of upper limits
; on [NII]/Ha and/or [OIII]/Hb
    ngal = n_elements(iratios1)
    iratios = struct_addtags(iratios1,replicate($
      {nii_ha_limit: -999.0, oiii_hb_limit: -999.0},ngal))
    nii = where((mz_ispec.h_alpha[1] gt 0.0) and (mz_ispec.nii_6584[1] eq -1.0),nnii)
    oiii = where((mz_ispec.oiii_5007[1] eq -1.0),noiii)
    iratios[nii].nii_ha_limit = alog10(mz_ispec[nii].nii_6584_limit/mz_ispec[nii].h_alpha[0])
    iratios[oiii].oiii_hb_limit = alog10(mz_ispec[oiii].oiii_5007_limit/mz_ispec[oiii].h_beta[0])

    case1 = where((strtrim(iratios.final_class,2) eq 'UNKNOWN') and $ ; [NII] limit, [OIII] limit
      (iratios.nii_ha_limit gt -900.0) and (iratios.oiii_hb_limit gt -900.0),ncase1)
    case2 = where((strtrim(iratios.final_class,2) eq 'UNKNOWN') and $ ; [NII] limit, [OIII] good
      (iratios.nii_ha_limit gt -900.0) and (iratios.oiii_hb gt -900.0),ncase2)
    case3 = where((strtrim(iratios.final_class,2) eq 'UNKNOWN') and $ ; [NII] good, [OIII] limit
      (iratios.nii_ha gt -900.0) and (iratios.oiii_hb_limit gt -900.0),ncase3)
    splog, 'niilimit_oiiilimit, niilimit_oiiigood, niigood_oiiilimit', ncase1, ncase2, ncase3

    if keyword_set(doplot) then begin
       djs_oplot, iratios[case1].nii_ha_limit, iratios[case1].oiii_hb_limit, psym=6, sym=0.2, color='orange'
       djs_oplot, iratios[case2].nii_ha_limit, iratios[case2].oiii_hb, psym=6, sym=0.2, color='cyan'
       djs_oplot, iratios[case3].nii_ha, iratios[case3].oiii_hb_limit, psym=6, sym=0.2, color='green'
    endif

; classify; note that the cut on max(kauffmann.x_nii) is necessary
; because in some cases the hyperbolic form of the demarcation line
; leads to some objects being misclassified
    kauffmann = kewley_bpt_lines(/kauffmann)

    data = kewley_bpt_lines(ratio_nii=iratios[case1].nii_ha_limit,/nii,/kauffmann)
    case1_sf = where((iratios[case1].oiii_hb_limit lt data.y_nii) and $
      (iratios[case1].nii_ha_limit lt max(kauffmann.x_nii)),ncase1_sf)
    iratios[case1[case1_sf]].final_class = 'BPT_SF_NIILIMIT_OIIILIMIT'

    data = kewley_bpt_lines(ratio_nii=iratios[case2].nii_ha_limit,/nii,/kauffmann)
    case2_sf = where((iratios[case2].oiii_hb lt data.y_nii) and $
      (iratios[case2].nii_ha_limit lt max(kauffmann.x_nii)),ncase2_sf)
    iratios[case2[case2_sf]].final_class = 'BPT_SF_NIILIMIT_OIIIGOOD'

    data = kewley_bpt_lines(ratio_nii=iratios[case3].nii_ha,/nii,/kauffmann)
    case3_sf = where((iratios[case3].oiii_hb_limit lt data.y_nii) and $
      (iratios[case3].nii_ha lt max(kauffmann.x_nii)),ncase3_sf) ; <-- need this cut on max()!
    iratios[case3[case3_sf]].final_class = 'BPT_SF_NIIGOOD_OIIILIMIT'

; debugging plot       
    check1 = where((strtrim(iratios[case1].final_class,2) eq 'UNKNOWN'))
    check2 = where((strtrim(iratios[case2].final_class,2) eq 'UNKNOWN'))
    check3 = where((strtrim(iratios[case3].final_class,2) eq 'UNKNOWN'))
;   djs_oplot, iratios[case1[check1]].nii_ha_limit, iratios[case1[check1]].oiii_hb_limit, psym=7, color='yellow'
;   djs_oplot, iratios[case2[check2]].nii_ha_limit, iratios[case2[check2]].oiii_hb, psym=7, color='yellow'
;   djs_oplot, iratios[case3[check3]].nii_ha, iratios[case3[check3]].oiii_hb_limit, psym=7, color='yellow'

; finally, we can classify CASE3 objects if they have
; log([NII]/Ha)>-0.3
    case3_agn = where((strtrim(iratios[case3].final_class,2) eq 'UNKNOWN') and $
      (iratios[case3].nii_ha gt -0.3),ncase3_agn)
    iratios[case3[case3_agn]].final_class = 'BPT_AGN_NIIGOOD_OIIILIMIT'
    if keyword_set(doplot) then begin
       djs_oplot, -0.3*[1,1], !y.crange, line=1
       djs_oplot, iratios[case3[case3_agn]].nii_ha, $
         iratios[case3[case3_agn]].oiii_hb_limit, psym=7, color='purple'
    endif

return, iratios
end

