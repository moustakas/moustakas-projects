pro write_sdss_abundances_sample, allsdssdust, allsdssnodust, allsdssancillary, write=write
; jm06apr16uofa - based on WRITE_SDSS_MZ_SAMPLEn

    outpath = atlas_path(/projects)+'abundances/'

    snrcut = 3.0
    ewhbcut = 5.0
    
    if (n_elements(allsdssdust) eq 0L) then allsdssdust = read_sdss(sdssnodust=allsdssnodust,ancillary=allsdssancillary)
    ngalaxy = n_elements(allsdssdust)

; define additional structures    

    ohdust = {$
      r23branch_pt05:              '?',$
      r23branch_pt05_ew:           '?',$
      r23branch_kk04:              '?',$
      r23branch_kk04_ew:           '?',$
      r23branch_m91:               '?',$
      r23branch_m91_ew:            '?',$

      zstrong_12oh_pt05:        -999.0,$
      zstrong_12oh_pt05_err:    -999.0,$
      zstrong_ew_12oh_pt05:     -999.0,$
      zstrong_ew_12oh_pt05_err: -999.0,$

      zstrong_12oh_m91:         -999.0,$
      zstrong_12oh_m91_err:     -999.0,$
      zstrong_ew_12oh_m91:      -999.0,$
      zstrong_ew_12oh_m91_err:  -999.0,$

      zstrong_12oh_kk04:        -999.0,$
      zstrong_12oh_kk04_err:    -999.0,$
      zstrong_ew_12oh_kk04:     -999.0,$
      zstrong_ew_12oh_kk04_err: -999.0}
    ohdust = replicate(ohdust,ngalaxy)
    ohnodust = ohdust
                    
; parent sample cuts

    parent = where((allsdssancillary.z_obj gt 0.033) and (allsdssancillary.sdss_g gt -900.0) and $
      (allsdssancillary.sdss_r gt -900.0) and (allsdssancillary.sdss_i gt -900.0),nparent)
;   nparent = n_elements(allsdssdust) & parent = lindgen(nparent)
;   parent = where((allsdssancillary.infiber_sdss_r gt 0.1) and (allsdssancillary.infiber_sdss_r le 1.0),nparent)
    splog, 'Parent sample: '+string(nparent,format='(I0)')+'/'+string(ngalaxy,format='(I0)')+' ('+$
      string(round(100.0*nparent/float(ngalaxy)),format='(I0)')+'%).'
    
;   mz = where((allsdssancillary.z_obj gt 0.033) and (allsdssancillary.sdss_g gt -900.0) and $
;     (allsdssancillary.sdss_r gt -900.0) and (allsdssancillary.sdss_i gt -900.0) and $
;     (allsdssdust.h_beta[0]/allsdssdust.h_beta[1] gt snrcut) and $
;     (allsdssdust.oii_3727[0]/allsdssdust.oii_3727[1] gt snrcut) and $
;     (allsdssdust.oiii_5007[0]/allsdssdust.oiii_5007[1] gt snrcut),nmz)
    mz = where((allsdssancillary.z_obj gt 0.033) and (allsdssancillary.sdss_g gt -900.0) and $
      (allsdssancillary.sdss_r gt -900.0) and (allsdssancillary.sdss_i gt -900.0) and $
      (allsdssdust.h_beta[0]/allsdssdust.h_beta[1] gt snrcut) and $
      (allsdssdust.h_alpha[0]/allsdssdust.h_alpha[1] gt snrcut) and $
      (allsdssdust.nii_6584[0]/allsdssdust.nii_6584[1] gt snrcut) and $
      (allsdssdust.oii_3727[0]/allsdssdust.oii_3727[1] gt snrcut) and $
      (allsdssdust.oiii_5007[0]/allsdssdust.oiii_5007[1] gt snrcut),nmz)

    hii = where(strtrim(allsdssdust[mz].bpt_pure_nii_class,2) eq 'HII',nhii)
    agn = where(strtrim(allsdssdust[mz].bpt_pure_nii_class,2) eq 'AGN',nagn)
    splog, 'MZ sample: '+string(nmz,format='(I0)')+'/'+string(nparent,format='(I0)')+' ('+$
      string(round(100.0*nmz/float(nparent)),format='(I0)')+'%).'
    splog, 'MZ sample (star-forming): '+string(nhii,format='(I0)')+'/'+string(nmz,format='(I0)')+' ('+$
      string(round(100.0*nhii/float(nmz)),format='(I0)')+'%).'
    splog, 'MZ sample (type 2 AGN): '+string(nagn,format='(I0)')+'/'+string(nmz,format='(I0)')+' ('+$
      string(round(100.0*nagn/float(nmz)),format='(I0)')+'%).'
    print

; figure out the appropriate R23 branch; the EW abundances are the
; same for ALLSDSSDUST and ALLSDSSNODUST; for now, assign the
; unassigned objects to the upper branch, but check back later for an
; LZ method

; PT05/lower branch    
    
    lo_pt05_dust = where((allsdssdust[mz].zstrong_niiha gt -900) and (allsdssdust[mz].zstrong_niiha le -1.0) and $
                         (allsdssdust[mz].zstrong_12oh_pt05_lower gt -900) and (allsdssdust[mz].zstrong_12oh_pt05_upper gt -900) and $
                         (allsdssdust[mz].zstrong_12oh_pt05_lower lt allsdssdust[mz].zstrong_12oh_pt05_upper),nlo_pt05_dust)
    if (nlo_pt05_dust ne 0L) then begin
       ohdust[mz[lo_pt05_dust]].r23branch_pt05   = 'L'
       ohdust[mz[lo_pt05_dust]].zstrong_12oh_pt05          = allsdssdust[mz[lo_pt05_dust]].zstrong_12oh_pt05_lower
       ohdust[mz[lo_pt05_dust]].zstrong_12oh_pt05_err      = allsdssdust[mz[lo_pt05_dust]].zstrong_12oh_pt05_lower_err
    endif

    lo_pt05_nodust = where((allsdssdust[mz].zstrong_niiha gt -900) and (allsdssdust[mz].zstrong_niiha le -1.0) and $
                           (allsdssnodust[mz].zstrong_12oh_pt05_lower gt -900) and (allsdssnodust[mz].zstrong_12oh_pt05_upper gt -900) and $
                           (allsdssnodust[mz].zstrong_12oh_pt05_lower lt allsdssnodust[mz].zstrong_12oh_pt05_upper),nlo_pt05_nodust)
    if (nlo_pt05_nodust ne 0L) then begin
       ohnodust[mz[lo_pt05_nodust]].r23branch_pt05           = 'L'
       ohnodust[mz[lo_pt05_nodust]].zstrong_12oh_pt05        = allsdssnodust[mz[lo_pt05_nodust]].zstrong_12oh_pt05_lower
       ohnodust[mz[lo_pt05_nodust]].zstrong_12oh_pt05_err    = allsdssnodust[mz[lo_pt05_nodust]].zstrong_12oh_pt05_lower_err
    endif

    lo_pt05_ew = where((allsdssdust[mz].zstrong_niiha gt -900) and (allsdssdust[mz].zstrong_niiha le -1.0) and $
                       (allsdssdust[mz].zstrong_ew_12oh_pt05_lower gt -900) and (allsdssdust[mz].zstrong_ew_12oh_pt05_upper gt -900) and $
                       (allsdssdust[mz].zstrong_ew_12oh_pt05_lower lt allsdssdust[mz].zstrong_ew_12oh_pt05_upper),nlo_pt05_ew)
    if (nlo_pt05_ew ne 0L) then begin
       ohdust[mz[lo_pt05_ew]].r23branch_pt05_ew            = 'L'
       ohdust[mz[lo_pt05_ew]].zstrong_ew_12oh_pt05         = allsdssdust[mz[lo_pt05_ew]].zstrong_ew_12oh_pt05_lower
       ohdust[mz[lo_pt05_ew]].zstrong_ew_12oh_pt05_err     = allsdssdust[mz[lo_pt05_ew]].zstrong_ew_12oh_pt05_lower_err
       ohnodust[mz[lo_pt05_ew]].r23branch_pt05_ew          = 'L'
       ohnodust[mz[lo_pt05_ew]].zstrong_ew_12oh_pt05       = allsdssdust[mz[lo_pt05_ew]].zstrong_ew_12oh_pt05_lower
       ohnodust[mz[lo_pt05_ew]].zstrong_ew_12oh_pt05_err   = allsdssdust[mz[lo_pt05_ew]].zstrong_ew_12oh_pt05_lower_err
    endif

; PT05/upper branch    

    up_pt05_dust = where((allsdssdust[mz].zstrong_niiha gt -900) and (allsdssdust[mz].zstrong_niiha gt -1.0) and $
                         (allsdssdust[mz].zstrong_12oh_pt05_upper gt -900) and (allsdssdust[mz].zstrong_12oh_pt05_lower gt -900) and $
                         (allsdssdust[mz].zstrong_12oh_pt05_upper gt allsdssdust[mz].zstrong_12oh_pt05_lower),nup_pt05_dust)
    if (nup_pt05_dust ne 0L) then begin
       ohdust[mz[up_pt05_dust]].r23branch_pt05             = 'U'
       ohdust[mz[up_pt05_dust]].zstrong_12oh_pt05          = allsdssdust[mz[up_pt05_dust]].zstrong_12oh_pt05_upper
       ohdust[mz[up_pt05_dust]].zstrong_12oh_pt05_err      = allsdssdust[mz[up_pt05_dust]].zstrong_12oh_pt05_upper_err
    endif

    up_pt05_nodust = where((allsdssdust[mz].zstrong_niiha gt -900) and (allsdssdust[mz].zstrong_niiha gt -1.0) and $
                           (allsdssnodust[mz].zstrong_12oh_pt05_upper gt -900) and (allsdssnodust[mz].zstrong_12oh_pt05_lower gt -900) and $
                           (allsdssnodust[mz].zstrong_12oh_pt05_upper gt allsdssnodust[mz].zstrong_12oh_pt05_lower),nup_pt05_nodust)
    if (nup_pt05_nodust ne 0L) then begin
       ohnodust[mz[up_pt05_nodust]].r23branch_pt05           = 'U'
       ohnodust[mz[up_pt05_nodust]].zstrong_12oh_pt05        = allsdssnodust[mz[up_pt05_nodust]].zstrong_12oh_pt05_upper
       ohnodust[mz[up_pt05_nodust]].zstrong_12oh_pt05_err    = allsdssnodust[mz[up_pt05_nodust]].zstrong_12oh_pt05_upper_err
    endif

    up_pt05_ew = where((allsdssdust[mz].zstrong_niiha gt -900) and (allsdssdust[mz].zstrong_niiha gt -1.0) and $
                       (allsdssdust[mz].zstrong_ew_12oh_pt05_lower gt -900) and (allsdssdust[mz].zstrong_ew_12oh_pt05_upper gt -900) and $
                       (allsdssdust[mz].zstrong_ew_12oh_pt05_upper gt allsdssdust[mz].zstrong_ew_12oh_pt05_lower),nup_pt05_ew)
    if (nup_pt05_ew ne 0L) then begin
       ohdust[mz[up_pt05_ew]].r23branch_pt05_ew            = 'U'
       ohdust[mz[up_pt05_ew]].zstrong_ew_12oh_pt05         = allsdssdust[mz[up_pt05_ew]].zstrong_ew_12oh_pt05_upper
       ohdust[mz[up_pt05_ew]].zstrong_ew_12oh_pt05_err     = allsdssdust[mz[up_pt05_ew]].zstrong_ew_12oh_pt05_upper_err
       ohnodust[mz[up_pt05_ew]].r23branch_pt05_ew          = 'U'
       ohnodust[mz[up_pt05_ew]].zstrong_ew_12oh_pt05       = allsdssdust[mz[up_pt05_ew]].zstrong_ew_12oh_pt05_upper
       ohnodust[mz[up_pt05_ew]].zstrong_ew_12oh_pt05_err   = allsdssdust[mz[up_pt05_ew]].zstrong_ew_12oh_pt05_upper_err
    endif
    
; PT05/ambigiuos    

    ambig_pt05_dust = where((allsdssdust[mz].zstrong_12oh_pt05_upper gt -900) and (allsdssdust[mz].zstrong_12oh_pt05_lower gt -900) and $
                            (allsdssdust[mz].zstrong_12oh_pt05_upper lt allsdssdust[mz].zstrong_12oh_pt05_lower) and $
                            ((allsdssdust[mz].zstrong_12oh_pt05_upper+allsdssdust[mz].zstrong_12oh_pt05_upper_err) gt $
                             (allsdssdust[mz].zstrong_12oh_pt05_lower-allsdssdust[mz].zstrong_12oh_pt05_lower_err)),nambig_pt05_dust)

;   niceprint, (allsdssdust[mz[ambig_pt05_dust]].zstrong_12oh_pt05_upper+allsdssdust[mz[ambig_pt05_dust]].zstrong_12oh_pt05_upper_err), $
;     (allsdssdust[mz[ambig_pt05_dust]].zstrong_12oh_pt05_lower-allsdssdust[mz[ambig_pt05_dust]].zstrong_12oh_pt05_lower_err), $
;     allsdssdust[mz[ambig_pt05_dust]].zstrong_12oh_pt05_upper_err, allsdssdust[mz[ambig_pt05_dust]].zstrong_12oh_pt05_lower_err
;   plot, allsdssdust[mz[ambig_pt05_dust]].zstrong_12oh_pt05_lower, allsdssdust[mz[ambig_pt05_dust]].zstrong_12oh_pt05_upper, $
;     ps=4, xsty=3, ysty=3, xrange=[7.5,9.0], yrange=[7.5,9.0] & oplot, !x.crange, !y.crange
;   ploterror, allsdssdust[mz[ambig_pt05_dust]].zstrong_12oh_pt05_lower, allsdssdust[mz[ambig_pt05_dust]].zstrong_12oh_pt05_upper, $
;     allsdssdust[mz[ambig_pt05_dust]].zstrong_12oh_pt05_lower_err, allsdssdust[mz[ambig_pt05_dust]].zstrong_12oh_pt05_upper_err, $
;     ps=4, xsty=3, ysty=3, xrange=[7.5,9.0], yrange=[7.5,9.0], errstyle=1 & oplot, !x.crange, !y.crange

;   ambig_pt05_dust = where((allsdssdust[mz].zstrong_12oh_pt05_upper gt -900) and (allsdssdust[mz].zstrong_12oh_pt05_lower gt -900) and $
;                           ((allsdssdust[mz].zstrong_12oh_pt05_upper lt allsdssdust[mz].zstrong_12oh_pt05_lower) or $
;                            (allsdssdust[mz].zstrong_12oh_pt05_lower gt allsdssdust[mz].zstrong_12oh_pt05_upper)) and $
;                           (abs(allsdssdust[mz].zstrong_12oh_pt05_upper-allsdssdust[mz].zstrong_12oh_pt05_lower) lt deltaoh),nambig_pt05_dust)

    if (nambig_pt05_dust ne 0L) then begin
       ohdust[mz[ambig_pt05_dust]].r23branch_pt05 = 'A'
       for i = 0L, nambig_pt05_dust-1L do begin
          oh     = [allsdssdust[mz[ambig_pt05_dust[i]]].zstrong_12oh_pt05_upper,$
                    allsdssdust[mz[ambig_pt05_dust[i]]].zstrong_12oh_pt05_lower]
          oh_err = [allsdssdust[mz[ambig_pt05_dust[i]]].zstrong_12oh_pt05_upper_err,$
                    allsdssdust[mz[ambig_pt05_dust[i]]].zstrong_12oh_pt05_lower_err]
          ohdust[mz[ambig_pt05_dust[i]]].zstrong_12oh_pt05     = total(oh/oh_err^2)/total(1.0/oh_err^2)
          ohdust[mz[ambig_pt05_dust[i]]].zstrong_12oh_pt05_err = 1.0/sqrt(total(1.0/oh_err^2))
       endfor
    endif

;   ploterror, ohdust[mz[ambig_pt05_dust]].zstrong_12oh_pt05, allsdssdust[mz[ambig_pt05_dust]].zstrong_12oh_pt05_lower, $
;     allsdssdust[mz[ambig_pt05_dust]].zstrong_12oh_pt05_lower_err, ps=4, xsty=3, ysty=3, xrange=[7.3,8.6], yrange=[7.3,8.6], $
;     errstyle=1
;   oploterror, ohdust[mz[ambig_pt05_dust]].zstrong_12oh_pt05, allsdssdust[mz[ambig_pt05_dust]].zstrong_12oh_pt05_upper, $
;     allsdssdust[mz[ambig_pt05_dust]].zstrong_12oh_pt05_upper_err, ps=4, errstyle=1, color=djs_icolor('red'), errcolor=djs_icolor('red')
;   oplot, !x.crange, !y.crange
    
    ambig_pt05_nodust = where((allsdssnodust[mz].zstrong_12oh_pt05_upper gt -900) and (allsdssnodust[mz].zstrong_12oh_pt05_lower gt -900) and $
                              (allsdssnodust[mz].zstrong_12oh_pt05_upper lt allsdssnodust[mz].zstrong_12oh_pt05_lower) and $
                              ((allsdssnodust[mz].zstrong_12oh_pt05_upper+allsdssnodust[mz].zstrong_12oh_pt05_upper_err) gt $
                               (allsdssnodust[mz].zstrong_12oh_pt05_lower-allsdssnodust[mz].zstrong_12oh_pt05_lower_err)),nambig_pt05_nodust)

    if (nambig_pt05_nodust ne 0L) then begin
       ohnodust[mz[ambig_pt05_nodust]].r23branch_pt05 = 'A'
       for i = 0L, nambig_pt05_nodust-1L do begin
          oh     = [allsdssnodust[mz[ambig_pt05_nodust[i]]].zstrong_12oh_pt05_upper,$
                    allsdssnodust[mz[ambig_pt05_nodust[i]]].zstrong_12oh_pt05_lower]
          oh_err = [allsdssnodust[mz[ambig_pt05_nodust[i]]].zstrong_12oh_pt05_upper_err,$
                    allsdssnodust[mz[ambig_pt05_nodust[i]]].zstrong_12oh_pt05_lower_err]
          ohnodust[mz[ambig_pt05_nodust[i]]].zstrong_12oh_pt05     = total(oh/oh_err^2)/total(1.0/oh_err^2)
          ohnodust[mz[ambig_pt05_nodust[i]]].zstrong_12oh_pt05_err = 1.0/sqrt(total(1.0/oh_err^2))
       endfor
    endif

    ambig_pt05_ew = where((allsdssdust[mz].zstrong_ew_12oh_pt05_upper gt -900) and (allsdssdust[mz].zstrong_ew_12oh_pt05_lower gt -900) and $
                          (allsdssdust[mz].zstrong_ew_12oh_pt05_upper lt allsdssdust[mz].zstrong_ew_12oh_pt05_lower) and $
                          ((allsdssdust[mz].zstrong_ew_12oh_pt05_upper+allsdssdust[mz].zstrong_ew_12oh_pt05_upper_err) gt $
                           (allsdssdust[mz].zstrong_ew_12oh_pt05_lower-allsdssdust[mz].zstrong_ew_12oh_pt05_lower_err)),nambig_pt05_ew)

    if (nambig_pt05_ew ne 0L) then begin
       ohdust[mz[ambig_pt05_ew]].r23branch_pt05_ew   = 'A'
       ohnodust[mz[ambig_pt05_ew]].r23branch_pt05_ew = 'A'
       for i = 0L, nambig_pt05_ew-1L do begin
          oh     = [allsdssdust[mz[ambig_pt05_ew[i]]].zstrong_ew_12oh_pt05_upper,$
                    allsdssdust[mz[ambig_pt05_ew[i]]].zstrong_ew_12oh_pt05_lower]
          oh_err = [allsdssdust[mz[ambig_pt05_ew[i]]].zstrong_ew_12oh_pt05_upper_err,$
                    allsdssdust[mz[ambig_pt05_ew[i]]].zstrong_ew_12oh_pt05_lower_err]
          ohdust[mz[ambig_pt05_ew[i]]].zstrong_ew_12oh_pt05       = total(oh/oh_err^2)/total(1.0/oh_err^2)
          ohdust[mz[ambig_pt05_ew[i]]].zstrong_ew_12oh_pt05_err   = 1.0/sqrt(total(1.0/oh_err^2))
          ohnodust[mz[ambig_pt05_ew[i]]].zstrong_ew_12oh_pt05     = total(oh/oh_err^2)/total(1.0/oh_err^2)
          ohnodust[mz[ambig_pt05_ew[i]]].zstrong_ew_12oh_pt05_err = 1.0/sqrt(total(1.0/oh_err^2))
       endfor
    endif

; for now, assign everything else to the PT05 upper branch

    need_pt05_dust = where((allsdssdust[mz].zstrong_12oh_pt05_upper gt -900.0) and (allsdssdust[mz].zstrong_12oh_pt05_lower gt -900.0) and $
                           (allsdssdust[mz].zstrong_12oh_pt05_upper gt allsdssdust[mz].zstrong_12oh_pt05_lower) and $
                           (ohdust[mz].zstrong_12oh_pt05 lt -900.0),nneed_pt05_dust)
    if (nneed_pt05_dust ne 0L) then begin
       ohdust[mz[need_pt05_dust]].r23branch_pt05             = 'U'
       ohdust[mz[need_pt05_dust]].zstrong_12oh_pt05          = allsdssdust[mz[need_pt05_dust]].zstrong_12oh_pt05_upper
       ohdust[mz[need_pt05_dust]].zstrong_12oh_pt05_err      = allsdssdust[mz[need_pt05_dust]].zstrong_12oh_pt05_upper_err
    endif

    need_pt05_nodust = where((allsdssnodust[mz].zstrong_12oh_pt05_upper gt -900.0) and (ohnodust[mz].zstrong_12oh_pt05 lt -900.0) and $
                             (allsdssnodust[mz].zstrong_12oh_pt05_upper gt allsdssnodust[mz].zstrong_12oh_pt05_lower) and $
                             (ohnodust[mz].zstrong_12oh_pt05 lt -900.0),nneed_pt05_nodust)
    if (nneed_pt05_nodust ne 0L) then begin
       ohnodust[mz[need_pt05_nodust]].r23branch_pt05             = 'U'
       ohnodust[mz[need_pt05_nodust]].zstrong_12oh_pt05          = allsdssnodust[mz[need_pt05_nodust]].zstrong_12oh_pt05_upper
       ohnodust[mz[need_pt05_nodust]].zstrong_12oh_pt05_err      = allsdssnodust[mz[need_pt05_nodust]].zstrong_12oh_pt05_upper_err
    endif

    need_pt05_ew_dust = where((allsdssdust[mz].zstrong_ew_12oh_pt05_upper gt -900.0) and (allsdssdust[mz].zstrong_ew_12oh_pt05_lower gt -900.0) and $
                              (allsdssdust[mz].zstrong_ew_12oh_pt05_upper gt allsdssdust[mz].zstrong_ew_12oh_pt05_lower) and $
                              (ohdust[mz].zstrong_ew_12oh_pt05 lt -900.0),nneed_pt05_ew_dust)
    if (nneed_pt05_ew_dust ne 0L) then begin
       ohdust[mz[need_pt05_ew_dust]].r23branch_pt05_ew        = 'U'
       ohdust[mz[need_pt05_ew_dust]].zstrong_ew_12oh_pt05     = allsdssdust[mz[need_pt05_ew_dust]].zstrong_ew_12oh_pt05_upper
       ohdust[mz[need_pt05_ew_dust]].zstrong_ew_12oh_pt05_err = allsdssdust[mz[need_pt05_ew_dust]].zstrong_ew_12oh_pt05_upper_err
    endif

    need_pt05_ew_nodust = where((allsdssdust[mz].zstrong_ew_12oh_pt05_upper gt -900.0) and (allsdssdust[mz].zstrong_ew_12oh_pt05_lower gt -900.0) and $
                                (allsdssdust[mz].zstrong_ew_12oh_pt05_upper gt allsdssdust[mz].zstrong_ew_12oh_pt05_lower) and $
                                (ohnodust[mz].zstrong_ew_12oh_pt05 lt -900.0),nneed_pt05_ew_nodust)
    if (nneed_pt05_ew_nodust ne 0L) then begin
       ohnodust[mz[need_pt05_ew_nodust]].r23branch_pt05_ew        = 'U'
       ohnodust[mz[need_pt05_ew_nodust]].zstrong_ew_12oh_pt05     = allsdssdust[mz[need_pt05_ew_nodust]].zstrong_ew_12oh_pt05_upper
       ohnodust[mz[need_pt05_ew_nodust]].zstrong_ew_12oh_pt05_err = allsdssdust[mz[need_pt05_ew_nodust]].zstrong_ew_12oh_pt05_upper_err
    endif

;   w = where((ohdust.zstrong_ew_12oh_pt05 gt -900) and (ohdust.r23branch_pt05_ew eq '?'),nw)
;   if (nw ne 0L) then stop
    
;   good = where((ohdust.zstrong_12oh_pt05 gt -900.0) and (ohdust.r23branch_pt05 eq '?'),ngood)
;   struct_print, struct_trimtags(ohdust[good],select='*PT05*')
;   if (ngood ne 0L) then stop

;   good = where(ohdust.zstrong_12oh_pt05 gt -900.0)
;   plot, allsdssdust[good].zstrong_r23, ohdust[good].zstrong_12oh_pt05, ps=4, xsty=3, ysty=3
;   struct_print, struct_trimtags(ohdust[good],select='*PT05*')

; M91/lower branch    
    
    lo_m91_dust = where((allsdssdust[mz].zstrong_niiha gt -900) and (allsdssdust[mz].zstrong_niiha le -1.0) and $
                         (allsdssdust[mz].zstrong_12oh_m91_lower gt -900) and (allsdssdust[mz].zstrong_12oh_m91_upper gt -900) and $
                         (allsdssdust[mz].zstrong_12oh_m91_lower lt allsdssdust[mz].zstrong_12oh_m91_upper),nlo_m91_dust)
    if (nlo_m91_dust ne 0L) then begin
       ohdust[mz[lo_m91_dust]].r23branch_m91   = 'L'
       ohdust[mz[lo_m91_dust]].zstrong_12oh_m91          = allsdssdust[mz[lo_m91_dust]].zstrong_12oh_m91_lower
       ohdust[mz[lo_m91_dust]].zstrong_12oh_m91_err      = allsdssdust[mz[lo_m91_dust]].zstrong_12oh_m91_lower_err
    endif

    lo_m91_nodust = where((allsdssdust[mz].zstrong_niiha gt -900) and (allsdssdust[mz].zstrong_niiha le -1.0) and $
                           (allsdssnodust[mz].zstrong_12oh_m91_lower gt -900) and (allsdssnodust[mz].zstrong_12oh_m91_upper gt -900) and $
                           (allsdssnodust[mz].zstrong_12oh_m91_lower lt allsdssnodust[mz].zstrong_12oh_m91_upper),nlo_m91_nodust)
    if (nlo_m91_nodust ne 0L) then begin
       ohnodust[mz[lo_m91_nodust]].r23branch_m91           = 'L'
       ohnodust[mz[lo_m91_nodust]].zstrong_12oh_m91        = allsdssnodust[mz[lo_m91_nodust]].zstrong_12oh_m91_lower
       ohnodust[mz[lo_m91_nodust]].zstrong_12oh_m91_err    = allsdssnodust[mz[lo_m91_nodust]].zstrong_12oh_m91_lower_err
    endif

    lo_m91_ew = where((allsdssdust[mz].zstrong_niiha gt -900) and (allsdssdust[mz].zstrong_niiha le -1.0) and $
                       (allsdssdust[mz].zstrong_ew_12oh_m91_lower gt -900) and (allsdssdust[mz].zstrong_ew_12oh_m91_upper gt -900) and $
                       (allsdssdust[mz].zstrong_ew_12oh_m91_lower lt allsdssdust[mz].zstrong_ew_12oh_m91_upper),nlo_m91_ew)
    if (nlo_m91_ew ne 0L) then begin
       ohdust[mz[lo_m91_ew]].r23branch_m91_ew            = 'L'
       ohdust[mz[lo_m91_ew]].zstrong_ew_12oh_m91         = allsdssdust[mz[lo_m91_ew]].zstrong_ew_12oh_m91_lower
       ohdust[mz[lo_m91_ew]].zstrong_ew_12oh_m91_err     = allsdssdust[mz[lo_m91_ew]].zstrong_ew_12oh_m91_lower_err
       ohnodust[mz[lo_m91_ew]].r23branch_m91_ew          = 'L'
       ohnodust[mz[lo_m91_ew]].zstrong_ew_12oh_m91       = allsdssdust[mz[lo_m91_ew]].zstrong_ew_12oh_m91_lower
       ohnodust[mz[lo_m91_ew]].zstrong_ew_12oh_m91_err   = allsdssdust[mz[lo_m91_ew]].zstrong_ew_12oh_m91_lower_err
    endif

; M91/upper branch    

    up_m91_dust = where((allsdssdust[mz].zstrong_niiha gt -900) and (allsdssdust[mz].zstrong_niiha gt -1.0) and $
                         (allsdssdust[mz].zstrong_12oh_m91_upper gt -900) and (allsdssdust[mz].zstrong_12oh_m91_lower gt -900) and $
                         (allsdssdust[mz].zstrong_12oh_m91_upper gt allsdssdust[mz].zstrong_12oh_m91_lower),nup_m91_dust)
    if (nup_m91_dust ne 0L) then begin
       ohdust[mz[up_m91_dust]].r23branch_m91             = 'U'
       ohdust[mz[up_m91_dust]].zstrong_12oh_m91          = allsdssdust[mz[up_m91_dust]].zstrong_12oh_m91_upper
       ohdust[mz[up_m91_dust]].zstrong_12oh_m91_err      = allsdssdust[mz[up_m91_dust]].zstrong_12oh_m91_upper_err
    endif

    up_m91_nodust = where((allsdssdust[mz].zstrong_niiha gt -900) and (allsdssdust[mz].zstrong_niiha gt -1.0) and $
                           (allsdssnodust[mz].zstrong_12oh_m91_upper gt -900) and (allsdssnodust[mz].zstrong_12oh_m91_lower gt -900) and $
                           (allsdssnodust[mz].zstrong_12oh_m91_upper gt allsdssnodust[mz].zstrong_12oh_m91_lower),nup_m91_nodust)
    if (nup_m91_nodust ne 0L) then begin
       ohnodust[mz[up_m91_nodust]].r23branch_m91           = 'U'
       ohnodust[mz[up_m91_nodust]].zstrong_12oh_m91        = allsdssnodust[mz[up_m91_nodust]].zstrong_12oh_m91_upper
       ohnodust[mz[up_m91_nodust]].zstrong_12oh_m91_err    = allsdssnodust[mz[up_m91_nodust]].zstrong_12oh_m91_upper_err
    endif

    up_m91_ew = where((allsdssdust[mz].zstrong_niiha gt -900) and (allsdssdust[mz].zstrong_niiha gt -1.0) and $
                       (allsdssdust[mz].zstrong_ew_12oh_m91_lower gt -900) and (allsdssdust[mz].zstrong_ew_12oh_m91_upper gt -900) and $
                       (allsdssdust[mz].zstrong_ew_12oh_m91_upper gt allsdssdust[mz].zstrong_ew_12oh_m91_lower),nup_m91_ew)
    if (nup_m91_ew ne 0L) then begin
       ohdust[mz[up_m91_ew]].r23branch_m91_ew            = 'U'
       ohdust[mz[up_m91_ew]].zstrong_ew_12oh_m91         = allsdssdust[mz[up_m91_ew]].zstrong_ew_12oh_m91_upper
       ohdust[mz[up_m91_ew]].zstrong_ew_12oh_m91_err     = allsdssdust[mz[up_m91_ew]].zstrong_ew_12oh_m91_upper_err
       ohnodust[mz[up_m91_ew]].r23branch_m91_ew          = 'U'
       ohnodust[mz[up_m91_ew]].zstrong_ew_12oh_m91       = allsdssdust[mz[up_m91_ew]].zstrong_ew_12oh_m91_upper
       ohnodust[mz[up_m91_ew]].zstrong_ew_12oh_m91_err   = allsdssdust[mz[up_m91_ew]].zstrong_ew_12oh_m91_upper_err
    endif
    
; M91/ambigiuos    

    ambig_m91_dust = where((allsdssdust[mz].zstrong_12oh_m91_upper gt -900) and (allsdssdust[mz].zstrong_12oh_m91_lower gt -900) and $
                            (allsdssdust[mz].zstrong_12oh_m91_upper lt allsdssdust[mz].zstrong_12oh_m91_lower) and $
                            ((allsdssdust[mz].zstrong_12oh_m91_upper+allsdssdust[mz].zstrong_12oh_m91_upper_err) gt $
                             (allsdssdust[mz].zstrong_12oh_m91_lower-allsdssdust[mz].zstrong_12oh_m91_lower_err)),nambig_m91_dust)

    if (nambig_m91_dust ne 0L) then begin
       ohdust[mz[ambig_m91_dust]].r23branch_m91 = 'A'
       for i = 0L, nambig_m91_dust-1L do begin
          oh     = [allsdssdust[mz[ambig_m91_dust[i]]].zstrong_12oh_m91_upper,$
                    allsdssdust[mz[ambig_m91_dust[i]]].zstrong_12oh_m91_lower]
          oh_err = [allsdssdust[mz[ambig_m91_dust[i]]].zstrong_12oh_m91_upper_err,$
                    allsdssdust[mz[ambig_m91_dust[i]]].zstrong_12oh_m91_lower_err]
          ohdust[mz[ambig_m91_dust[i]]].zstrong_12oh_m91     = total(oh/oh_err^2)/total(1.0/oh_err^2)
          ohdust[mz[ambig_m91_dust[i]]].zstrong_12oh_m91_err = 1.0/sqrt(total(1.0/oh_err^2))
       endfor
    endif

    ambig_m91_nodust = where((allsdssnodust[mz].zstrong_12oh_m91_upper gt -900) and (allsdssnodust[mz].zstrong_12oh_m91_lower gt -900) and $
                             (allsdssnodust[mz].zstrong_12oh_m91_upper lt allsdssnodust[mz].zstrong_12oh_m91_lower) and $
                             ((allsdssnodust[mz].zstrong_12oh_m91_upper+allsdssnodust[mz].zstrong_12oh_m91_upper_err) gt $
                              (allsdssnodust[mz].zstrong_12oh_m91_lower-allsdssnodust[mz].zstrong_12oh_m91_lower_err)),nambig_m91_nodust)

    if (nambig_m91_nodust ne 0L) then begin
       ohnodust[mz[ambig_m91_nodust]].r23branch_m91 = 'A'
       for i = 0L, nambig_m91_nodust-1L do begin
          oh     = [allsdssnodust[mz[ambig_m91_nodust[i]]].zstrong_12oh_m91_upper,$
                    allsdssnodust[mz[ambig_m91_nodust[i]]].zstrong_12oh_m91_lower]
          oh_err = [allsdssnodust[mz[ambig_m91_nodust[i]]].zstrong_12oh_m91_upper_err,$
                    allsdssnodust[mz[ambig_m91_nodust[i]]].zstrong_12oh_m91_lower_err]
          ohnodust[mz[ambig_m91_nodust[i]]].zstrong_12oh_m91     = total(oh/oh_err^2)/total(1.0/oh_err^2)
          ohnodust[mz[ambig_m91_nodust[i]]].zstrong_12oh_m91_err = 1.0/sqrt(total(1.0/oh_err^2))
       endfor
    endif

    ambig_m91_ew = where((allsdssdust[mz].zstrong_ew_12oh_m91_upper gt -900) and (allsdssdust[mz].zstrong_ew_12oh_m91_lower gt -900) and $
                          (allsdssdust[mz].zstrong_ew_12oh_m91_upper lt allsdssdust[mz].zstrong_ew_12oh_m91_lower) and $
                          ((allsdssdust[mz].zstrong_ew_12oh_m91_upper+allsdssdust[mz].zstrong_ew_12oh_m91_upper_err) gt $
                           (allsdssdust[mz].zstrong_ew_12oh_m91_lower-allsdssdust[mz].zstrong_ew_12oh_m91_lower_err)),nambig_m91_ew)

    if (nambig_m91_ew ne 0L) then begin
       ohdust[mz[ambig_m91_ew]].r23branch_m91_ew   = 'A'
       ohnodust[mz[ambig_m91_ew]].r23branch_m91_ew = 'A'
       for i = 0L, nambig_m91_ew-1L do begin
          oh     = [allsdssdust[mz[ambig_m91_ew[i]]].zstrong_ew_12oh_m91_upper,$
                    allsdssdust[mz[ambig_m91_ew[i]]].zstrong_ew_12oh_m91_lower]
          oh_err = [allsdssdust[mz[ambig_m91_ew[i]]].zstrong_ew_12oh_m91_upper_err,$
                    allsdssdust[mz[ambig_m91_ew[i]]].zstrong_ew_12oh_m91_lower_err]
          ohdust[mz[ambig_m91_ew[i]]].zstrong_ew_12oh_m91       = total(oh/oh_err^2)/total(1.0/oh_err^2)
          ohdust[mz[ambig_m91_ew[i]]].zstrong_ew_12oh_m91_err   = 1.0/sqrt(total(1.0/oh_err^2))
          ohnodust[mz[ambig_m91_ew[i]]].zstrong_ew_12oh_m91     = total(oh/oh_err^2)/total(1.0/oh_err^2)
          ohnodust[mz[ambig_m91_ew[i]]].zstrong_ew_12oh_m91_err = 1.0/sqrt(total(1.0/oh_err^2))
       endfor
    endif

; for now, assign everything else to the M91 upper branch

    need_m91_dust = where((allsdssdust[mz].zstrong_12oh_m91_upper gt -900.0) and (allsdssdust[mz].zstrong_12oh_m91_lower gt -900.0) and $
                           (allsdssdust[mz].zstrong_12oh_m91_upper gt allsdssdust[mz].zstrong_12oh_m91_lower) and $
                           (ohdust[mz].zstrong_12oh_m91 lt -900.0),nneed_m91_dust)
    if (nneed_m91_dust ne 0L) then begin
       ohdust[mz[need_m91_dust]].r23branch_m91             = 'U'
       ohdust[mz[need_m91_dust]].zstrong_12oh_m91          = allsdssdust[mz[need_m91_dust]].zstrong_12oh_m91_upper
       ohdust[mz[need_m91_dust]].zstrong_12oh_m91_err      = allsdssdust[mz[need_m91_dust]].zstrong_12oh_m91_upper_err
    endif

    need_m91_nodust = where((allsdssnodust[mz].zstrong_12oh_m91_upper gt -900.0) and (ohnodust[mz].zstrong_12oh_m91 lt -900.0) and $
                             (allsdssnodust[mz].zstrong_12oh_m91_upper gt allsdssnodust[mz].zstrong_12oh_m91_lower) and $
                             (ohnodust[mz].zstrong_12oh_m91 lt -900.0),nneed_m91_nodust)
    if (nneed_m91_nodust ne 0L) then begin
       ohnodust[mz[need_m91_nodust]].r23branch_m91             = 'U'
       ohnodust[mz[need_m91_nodust]].zstrong_12oh_m91          = allsdssnodust[mz[need_m91_nodust]].zstrong_12oh_m91_upper
       ohnodust[mz[need_m91_nodust]].zstrong_12oh_m91_err      = allsdssnodust[mz[need_m91_nodust]].zstrong_12oh_m91_upper_err
    endif

    need_m91_ew_dust = where((allsdssdust[mz].zstrong_ew_12oh_m91_upper gt -900.0) and (allsdssdust[mz].zstrong_ew_12oh_m91_lower gt -900.0) and $
                              (allsdssdust[mz].zstrong_ew_12oh_m91_upper gt allsdssdust[mz].zstrong_ew_12oh_m91_lower) and $
                              (ohdust[mz].zstrong_ew_12oh_m91 lt -900.0),nneed_m91_ew_dust)
    if (nneed_m91_ew_dust ne 0L) then begin
       ohdust[mz[need_m91_ew_dust]].r23branch_m91_ew        = 'U'
       ohdust[mz[need_m91_ew_dust]].zstrong_ew_12oh_m91     = allsdssdust[mz[need_m91_ew_dust]].zstrong_ew_12oh_m91_upper
       ohdust[mz[need_m91_ew_dust]].zstrong_ew_12oh_m91_err = allsdssdust[mz[need_m91_ew_dust]].zstrong_ew_12oh_m91_upper_err
    endif

    need_m91_ew_nodust = where((allsdssdust[mz].zstrong_ew_12oh_m91_upper gt -900.0) and (allsdssdust[mz].zstrong_ew_12oh_m91_lower gt -900.0) and $
                                (allsdssdust[mz].zstrong_ew_12oh_m91_upper gt allsdssdust[mz].zstrong_ew_12oh_m91_lower) and $
                                (ohnodust[mz].zstrong_ew_12oh_m91 lt -900.0),nneed_m91_ew_nodust)
    if (nneed_m91_ew_nodust ne 0L) then begin
       ohnodust[mz[need_m91_ew_nodust]].r23branch_m91_ew        = 'U'
       ohnodust[mz[need_m91_ew_nodust]].zstrong_ew_12oh_m91     = allsdssdust[mz[need_m91_ew_nodust]].zstrong_ew_12oh_m91_upper
       ohnodust[mz[need_m91_ew_nodust]].zstrong_ew_12oh_m91_err = allsdssdust[mz[need_m91_ew_nodust]].zstrong_ew_12oh_m91_upper_err
    endif

;   w = where((ohdust.zstrong_ew_12oh_m91 gt -900) and (ohdust.r23branch_m91_ew eq '?'),nw)
;   if (nw ne 0L) then stop
    
;   good = where((ohdust.zstrong_12oh_m91 gt -900.0) and (ohdust.r23branch_m91 eq '?'),ngood)
;   struct_print, struct_trimtags(ohdust[good],select='*M91*')
;   if (ngood ne 0L) then stop

;   good = where(ohdust.zstrong_12oh_m91 gt -900.0)
;   plot, allsdssdust[good].zstrong_r23, ohdust[good].zstrong_12oh_m91, ps=4, xsty=3, ysty=3
;   struct_print, struct_trimtags(ohdust[good],select='*M91*')

; KK04/lower branch    
    
    lo_kk04_dust = where((allsdssdust[mz].zstrong_niiha gt -900) and (allsdssdust[mz].zstrong_niiha le -1.0) and $
                         (allsdssdust[mz].zstrong_12oh_kk04_r23_lower gt -900) and (allsdssdust[mz].zstrong_12oh_kk04_r23_upper gt -900) and $
                         (allsdssdust[mz].zstrong_12oh_kk04_r23_lower lt allsdssdust[mz].zstrong_12oh_kk04_r23_upper),nlo_kk04_dust)
    if (nlo_kk04_dust ne 0L) then begin
       ohdust[mz[lo_kk04_dust]].r23branch_kk04             = 'L'
       ohdust[mz[lo_kk04_dust]].zstrong_12oh_kk04          = allsdssdust[mz[lo_kk04_dust]].zstrong_12oh_kk04_r23_lower
       ohdust[mz[lo_kk04_dust]].zstrong_12oh_kk04_err      = allsdssdust[mz[lo_kk04_dust]].zstrong_12oh_kk04_r23_lower_err
    endif

    lo_kk04_nodust = where((allsdssdust[mz].zstrong_niiha gt -900) and (allsdssdust[mz].zstrong_niiha le -1.0) and $
                           (allsdssnodust[mz].zstrong_12oh_kk04_r23_lower gt -900) and (allsdssnodust[mz].zstrong_12oh_kk04_r23_upper gt -900) and $
                           (allsdssnodust[mz].zstrong_12oh_kk04_r23_lower lt allsdssnodust[mz].zstrong_12oh_kk04_r23_upper),nlo_kk04_nodust)
    if (nlo_kk04_nodust ne 0L) then begin
       ohnodust[mz[lo_kk04_nodust]].r23branch_kk04           = 'L'
       ohnodust[mz[lo_kk04_nodust]].zstrong_12oh_kk04        = allsdssnodust[mz[lo_kk04_nodust]].zstrong_12oh_kk04_r23_lower
       ohnodust[mz[lo_kk04_nodust]].zstrong_12oh_kk04_err    = allsdssnodust[mz[lo_kk04_nodust]].zstrong_12oh_kk04_r23_lower_err
    endif

    lo_kk04_ew = where((allsdssdust[mz].zstrong_niiha gt -900) and (allsdssdust[mz].zstrong_niiha le -1.0) and $
                       (allsdssdust[mz].zstrong_ew_12oh_kk04_r23_lower gt -900) and (allsdssdust[mz].zstrong_ew_12oh_kk04_r23_upper gt -900) and $
                       (allsdssdust[mz].zstrong_ew_12oh_kk04_r23_lower lt allsdssdust[mz].zstrong_ew_12oh_kk04_r23_upper),nlo_kk04_ew)
    if (nlo_kk04_ew ne 0L) then begin
       ohdust[mz[lo_kk04_ew]].r23branch_kk04_ew            = 'L'
       ohdust[mz[lo_kk04_ew]].zstrong_ew_12oh_kk04         = allsdssdust[mz[lo_kk04_ew]].zstrong_ew_12oh_kk04_r23_lower
       ohdust[mz[lo_kk04_ew]].zstrong_ew_12oh_kk04_err     = allsdssdust[mz[lo_kk04_ew]].zstrong_ew_12oh_kk04_r23_lower_err
       ohnodust[mz[lo_kk04_ew]].r23branch_kk04_ew          = 'L'
       ohnodust[mz[lo_kk04_ew]].zstrong_ew_12oh_kk04       = allsdssdust[mz[lo_kk04_ew]].zstrong_ew_12oh_kk04_r23_lower
       ohnodust[mz[lo_kk04_ew]].zstrong_ew_12oh_kk04_err   = allsdssdust[mz[lo_kk04_ew]].zstrong_ew_12oh_kk04_r23_lower_err
    endif
    
; KK04/upper branch    

    up_kk04_dust = where((allsdssdust[mz].zstrong_niiha gt -900) and (allsdssdust[mz].zstrong_niiha gt -1.0) and $
                         (allsdssdust[mz].zstrong_12oh_kk04_r23_upper gt -900) and (allsdssdust[mz].zstrong_12oh_kk04_r23_lower gt -900) and $
                         (allsdssdust[mz].zstrong_12oh_kk04_r23_upper gt allsdssdust[mz].zstrong_12oh_kk04_r23_lower),nup_kk04_dust)
    if (nup_kk04_dust ne 0L) then begin
       ohdust[mz[up_kk04_dust]].r23branch_kk04        = 'U'
       ohdust[mz[up_kk04_dust]].zstrong_12oh_kk04     = allsdssdust[mz[up_kk04_dust]].zstrong_12oh_kk04_r23_upper
       ohdust[mz[up_kk04_dust]].zstrong_12oh_kk04_err = allsdssdust[mz[up_kk04_dust]].zstrong_12oh_kk04_r23_upper_err
    endif

    up_kk04_nodust = where((allsdssdust[mz].zstrong_niiha gt -900) and (allsdssdust[mz].zstrong_niiha gt -1.0) and $
                           (allsdssnodust[mz].zstrong_12oh_kk04_r23_upper gt -900) and (allsdssnodust[mz].zstrong_12oh_kk04_r23_lower gt -900) and $
                           (allsdssnodust[mz].zstrong_12oh_kk04_r23_upper gt allsdssnodust[mz].zstrong_12oh_kk04_r23_lower),nup_kk04_nodust)
    if (nup_kk04_nodust ne 0L) then begin
       ohnodust[mz[up_kk04_nodust]].r23branch_kk04        = 'U'
       ohnodust[mz[up_kk04_nodust]].zstrong_12oh_kk04     = allsdssnodust[mz[up_kk04_nodust]].zstrong_12oh_kk04_r23_upper
       ohnodust[mz[up_kk04_nodust]].zstrong_12oh_kk04_err = allsdssnodust[mz[up_kk04_nodust]].zstrong_12oh_kk04_r23_upper_err
    endif

    up_kk04_ew = where((allsdssdust[mz].zstrong_niiha gt -900) and (allsdssdust[mz].zstrong_niiha gt -1.0) and $
                       (allsdssdust[mz].zstrong_ew_12oh_kk04_r23_lower gt -900) and (allsdssdust[mz].zstrong_ew_12oh_kk04_r23_upper gt -900) and $
                       (allsdssdust[mz].zstrong_ew_12oh_kk04_r23_upper gt allsdssdust[mz].zstrong_ew_12oh_kk04_r23_lower),nup_kk04_ew)
    if (nup_kk04_ew ne 0L) then begin
       ohdust[mz[up_kk04_ew]].r23branch_kk04_ew          = 'U'
       ohdust[mz[up_kk04_ew]].zstrong_ew_12oh_kk04       = allsdssdust[mz[up_kk04_ew]].zstrong_ew_12oh_kk04_r23_upper
       ohdust[mz[up_kk04_ew]].zstrong_ew_12oh_kk04_err   = allsdssdust[mz[up_kk04_ew]].zstrong_ew_12oh_kk04_r23_upper_err
       ohnodust[mz[up_kk04_ew]].r23branch_kk04_ew        = 'U'
       ohnodust[mz[up_kk04_ew]].zstrong_ew_12oh_kk04     = allsdssdust[mz[up_kk04_ew]].zstrong_ew_12oh_kk04_r23_upper
       ohnodust[mz[up_kk04_ew]].zstrong_ew_12oh_kk04_err = allsdssdust[mz[up_kk04_ew]].zstrong_ew_12oh_kk04_r23_upper_err
    endif

; KK04/ambigiuos    

    ambig_kk04_dust = where((allsdssdust[mz].zstrong_12oh_kk04_r23_upper gt -900) and (allsdssdust[mz].zstrong_12oh_kk04_r23_lower gt -900) and $
                            (allsdssdust[mz].zstrong_12oh_kk04_r23_upper lt allsdssdust[mz].zstrong_12oh_kk04_r23_lower) and $
                            ((allsdssdust[mz].zstrong_12oh_kk04_r23_upper+allsdssdust[mz].zstrong_12oh_kk04_r23_upper_err) gt $
                             (allsdssdust[mz].zstrong_12oh_kk04_r23_lower-allsdssdust[mz].zstrong_12oh_kk04_r23_lower_err)),nambig_kk04_dust)

    if (nambig_kk04_dust ne 0L) then begin
       ohdust[mz[ambig_kk04_dust]].r23branch_kk04 = 'A'
       for i = 0L, nambig_kk04_dust-1L do begin
          oh     = [allsdssdust[mz[ambig_kk04_dust[i]]].zstrong_12oh_kk04_r23_upper,$
                    allsdssdust[mz[ambig_kk04_dust[i]]].zstrong_12oh_kk04_r23_lower]
          oh_err = [allsdssdust[mz[ambig_kk04_dust[i]]].zstrong_12oh_kk04_r23_upper_err,$
                    allsdssdust[mz[ambig_kk04_dust[i]]].zstrong_12oh_kk04_r23_lower_err]
          ohdust[mz[ambig_kk04_dust[i]]].zstrong_12oh_kk04     = total(oh/oh_err^2)/total(1.0/oh_err^2)
          ohdust[mz[ambig_kk04_dust[i]]].zstrong_12oh_kk04_err = 1.0/sqrt(total(1.0/oh_err^2))
       endfor
    endif

    ambig_kk04_nodust = where((allsdssnodust[mz].zstrong_12oh_kk04_r23_upper gt -900) and (allsdssnodust[mz].zstrong_12oh_kk04_r23_lower gt -900) and $
                            (allsdssnodust[mz].zstrong_12oh_kk04_r23_upper lt allsdssnodust[mz].zstrong_12oh_kk04_r23_lower) and $
                            ((allsdssnodust[mz].zstrong_12oh_kk04_r23_upper+allsdssnodust[mz].zstrong_12oh_kk04_r23_upper_err) gt $
                             (allsdssnodust[mz].zstrong_12oh_kk04_r23_lower-allsdssnodust[mz].zstrong_12oh_kk04_r23_lower_err)),nambig_kk04_nodust)

    if (nambig_kk04_nodust ne 0L) then begin
       ohnodust[mz[ambig_kk04_nodust]].r23branch_kk04 = 'A'
       for i = 0L, nambig_kk04_nodust-1L do begin
          oh     = [allsdssnodust[mz[ambig_kk04_nodust[i]]].zstrong_12oh_kk04_r23_upper,$
                    allsdssnodust[mz[ambig_kk04_nodust[i]]].zstrong_12oh_kk04_r23_lower]
          oh_err = [allsdssnodust[mz[ambig_kk04_nodust[i]]].zstrong_12oh_kk04_r23_upper_err,$
                    allsdssnodust[mz[ambig_kk04_nodust[i]]].zstrong_12oh_kk04_r23_lower_err]
          ohnodust[mz[ambig_kk04_nodust[i]]].zstrong_12oh_kk04     = total(oh/oh_err^2)/total(1.0/oh_err^2)
          ohnodust[mz[ambig_kk04_nodust[i]]].zstrong_12oh_kk04_err = 1.0/sqrt(total(1.0/oh_err^2))
       endfor
    endif

    ambig_kk04_ew = where((allsdssdust[mz].zstrong_ew_12oh_kk04_r23_upper gt -900) and (allsdssdust[mz].zstrong_ew_12oh_kk04_r23_lower gt -900) and $
                          (allsdssdust[mz].zstrong_ew_12oh_kk04_r23_upper lt allsdssdust[mz].zstrong_ew_12oh_kk04_r23_lower) and $
                          ((allsdssdust[mz].zstrong_ew_12oh_kk04_r23_upper+allsdssdust[mz].zstrong_ew_12oh_kk04_r23_upper_err) gt $
                           (allsdssdust[mz].zstrong_ew_12oh_kk04_r23_lower-allsdssdust[mz].zstrong_ew_12oh_kk04_r23_lower_err)),nambig_kk04_ew)

    if (nambig_kk04_ew ne 0L) then begin
       ohdust[mz[ambig_kk04_ew]].r23branch_kk04_ew   = 'A'
       ohnodust[mz[ambig_kk04_ew]].r23branch_kk04_ew = 'A'
       for i = 0L, nambig_kk04_ew-1L do begin
          oh     = [allsdssdust[mz[ambig_kk04_ew[i]]].zstrong_ew_12oh_kk04_r23_upper,$
                    allsdssdust[mz[ambig_kk04_ew[i]]].zstrong_ew_12oh_kk04_r23_lower]
          oh_err = [allsdssdust[mz[ambig_kk04_ew[i]]].zstrong_ew_12oh_kk04_r23_upper_err,$
                    allsdssdust[mz[ambig_kk04_ew[i]]].zstrong_ew_12oh_kk04_r23_lower_err]
          ohdust[mz[ambig_kk04_ew[i]]].zstrong_ew_12oh_kk04       = total(oh/oh_err^2)/total(1.0/oh_err^2)
          ohdust[mz[ambig_kk04_ew[i]]].zstrong_ew_12oh_kk04_err   = 1.0/sqrt(total(1.0/oh_err^2))
          ohnodust[mz[ambig_kk04_ew[i]]].zstrong_ew_12oh_kk04     = total(oh/oh_err^2)/total(1.0/oh_err^2)
          ohnodust[mz[ambig_kk04_ew[i]]].zstrong_ew_12oh_kk04_err = 1.0/sqrt(total(1.0/oh_err^2))
       endfor
    endif

; for now, assign everything else to the upper branch

    need_kk04_dust = where((allsdssdust[mz].zstrong_12oh_kk04_r23_upper gt -900.0) and (allsdssdust[mz].zstrong_12oh_kk04_r23_lower gt -900.0) and $
                           (allsdssdust[mz].zstrong_12oh_kk04_r23_upper gt allsdssdust[mz].zstrong_12oh_kk04_r23_lower) and $
                           (ohdust[mz].zstrong_12oh_kk04 lt -900.0),nneed_kk04_dust)
    if (nneed_kk04_dust ne 0L) then begin
       ohdust[mz[need_kk04_dust]].r23branch_kk04        = 'U'
       ohdust[mz[need_kk04_dust]].zstrong_12oh_kk04     = allsdssdust[mz[need_kk04_dust]].zstrong_12oh_kk04_r23_upper
       ohdust[mz[need_kk04_dust]].zstrong_12oh_kk04_err = allsdssdust[mz[need_kk04_dust]].zstrong_12oh_kk04_r23_upper_err
    endif

    need_kk04_nodust = where((allsdssnodust[mz].zstrong_12oh_kk04_r23_upper gt -900.0) and (ohnodust[mz].zstrong_12oh_kk04 lt -900.0) and $
                             (allsdssnodust[mz].zstrong_12oh_kk04_r23_upper gt allsdssnodust[mz].zstrong_12oh_kk04_r23_lower) and $
                             (ohnodust[mz].zstrong_12oh_kk04 lt -900.0),nneed_kk04_nodust)
    if (nneed_kk04_nodust ne 0L) then begin
       ohnodust[mz[need_kk04_nodust]].r23branch_kk04        = 'U'
       ohnodust[mz[need_kk04_nodust]].zstrong_12oh_kk04     = allsdssnodust[mz[need_kk04_nodust]].zstrong_12oh_kk04_r23_upper
       ohnodust[mz[need_kk04_nodust]].zstrong_12oh_kk04_err = allsdssnodust[mz[need_kk04_nodust]].zstrong_12oh_kk04_r23_upper_err
    endif

    need_kk04_ew_dust = where((allsdssdust[mz].zstrong_ew_12oh_kk04_r23_upper gt -900.0) and (allsdssdust[mz].zstrong_ew_12oh_kk04_r23_lower gt -900.0) and $
                              (allsdssdust[mz].zstrong_ew_12oh_kk04_r23_upper gt allsdssdust[mz].zstrong_ew_12oh_kk04_r23_lower) and $
                              (ohdust[mz].zstrong_ew_12oh_kk04 lt -900.0),nneed_kk04_ew_dust)
    if (nneed_kk04_ew_dust ne 0L) then begin
       ohdust[mz[need_kk04_ew_dust]].r23branch_kk04_ew        = 'U'
       ohdust[mz[need_kk04_ew_dust]].zstrong_ew_12oh_kk04     = allsdssdust[mz[need_kk04_ew_dust]].zstrong_ew_12oh_kk04_r23_upper
       ohdust[mz[need_kk04_ew_dust]].zstrong_ew_12oh_kk04_err = allsdssdust[mz[need_kk04_ew_dust]].zstrong_ew_12oh_kk04_r23_upper_err
    endif

    need_kk04_ew_nodust = where((allsdssdust[mz].zstrong_ew_12oh_kk04_r23_upper gt -900.0) and (allsdssdust[mz].zstrong_ew_12oh_kk04_r23_lower gt -900.0) and $
                              (allsdssdust[mz].zstrong_ew_12oh_kk04_r23_upper gt allsdssdust[mz].zstrong_ew_12oh_kk04_r23_lower) and $
                              (ohnodust[mz].zstrong_ew_12oh_kk04 lt -900.0),nneed_kk04_ew_nodust)
    if (nneed_kk04_ew_nodust ne 0L) then begin
       ohnodust[mz[need_kk04_ew_nodust]].r23branch_kk04_ew        = 'U'
       ohnodust[mz[need_kk04_ew_nodust]].zstrong_ew_12oh_kk04     = allsdssdust[mz[need_kk04_ew_nodust]].zstrong_ew_12oh_kk04_r23_upper
       ohnodust[mz[need_kk04_ew_nodust]].zstrong_ew_12oh_kk04_err = allsdssdust[mz[need_kk04_ew_nodust]].zstrong_ew_12oh_kk04_r23_upper_err
    endif

;   w = where((ohdust.zstrong_ew_12oh_kk04 gt -900) and (ohdust.r23branch_kk04_ew eq '?'),nw)
;   if (nw ne 0L) then stop
    
;   good = where(ohdust.zstrong_12oh_kk04 gt -900.0)
;   plot, allsdssdust[good].zstrong_r23, ohdust[good].zstrong_12oh_kk04, ps=4, xsty=3, ysty=3, syms=0.5
;   good = where((ohdust.zstrong_12oh_kk04 gt -900.0) and (ohdust.zstrong_12oh_pt05 gt -900.0))
;   plot, ohdust[good].zstrong_12oh_kk04, ohdust[good].zstrong_12oh_pt05, ps=4, xsty=3, ysty=3, syms=0.5
;   good = where((ohnodust.zstrong_12oh_kk04 gt -900.0) and (ohnodust.zstrong_12oh_pt05 gt -900.0))
;   plot, ohnodust[good].zstrong_12oh_kk04, ohnodust[good].zstrong_12oh_pt05, ps=4, xsty=3, ysty=3, syms=0.5
;   good = where((ohdust.zstrong_ew_12oh_pt05 gt -900.0) and (ohdust.zstrong_12oh_pt05 gt -900.0))
;   plot, ohdust[good].zstrong_ew_12oh_pt05, ohdust[good].zstrong_12oh_pt05, ps=4, xsty=3, ysty=3, syms=0.5
;   good = where((ohdust.zstrong_ew_12oh_pt05 gt -900.0) and (ohnodust.zstrong_12oh_pt05 gt -900.0))
;   plot, ohdust[good].zstrong_ew_12oh_pt05, ohnodust[good].zstrong_12oh_pt05, ps=4, xsty=3, ysty=3, syms=0.5
;   good = where((ohdust.zstrong_ew_12oh_kk04 gt -900.0) and (ohdust.zstrong_12oh_kk04 gt -900.0))
;   plot, ohdust[good].zstrong_ew_12oh_kk04, ohdust[good].zstrong_12oh_kk04, ps=4, xsty=3, ysty=3, syms=0.5
;   good = where((ohdust.zstrong_ew_12oh_kk04 gt -900.0) and (ohnodust.zstrong_12oh_kk04 gt -900.0))
;   plot, ohdust[good].zstrong_ew_12oh_kk04, ohnodust[good].zstrong_12oh_kk04, ps=4, xsty=3, ysty=3, syms=0.5

; print some statistics

    pt05_good = where((ohdust[mz].zstrong_12oh_pt05 gt -900.0),npt05_good)
    pt05_good_upper = where((ohdust[mz[pt05_good]].r23branch_pt05 eq 'U'),npt05_good_upper)
    pt05_good_lower = where((ohdust[mz[pt05_good]].r23branch_pt05 eq 'L'),npt05_good_lower)
    pt05_good_ambig = where((ohdust[mz[pt05_good]].r23branch_pt05 eq 'A'),npt05_good_ambig)
    splog, 'PT05 abundances: '+string(npt05_good,format='(I0)')+'/'+string(nmz,format='(I0)')+' ('+$
      string(round(100.0*npt05_good/nmz),format='(I0)')+'%).'
    splog, '   Upper: '+string(npt05_good_upper,format='(I0)')+'/'+string(npt05_good,format='(I0)')+' ('+$
      string(round(100.0*npt05_good_upper/npt05_good),format='(I0)')+'%).'
    splog, '   Lower: '+string(npt05_good_lower,format='(I0)')+'/'+string(npt05_good,format='(I0)')+' ('+$
      string(round(100.0*npt05_good_lower/npt05_good),format='(I0)')+'%).'
    splog, '   Ambig: '+string(npt05_good_ambig,format='(I0)')+'/'+string(npt05_good,format='(I0)')+' ('+$
      string(round(100.0*npt05_good_ambig/npt05_good),format='(I0)')+'%).'
    print

    m91_good = where((ohdust[mz].zstrong_12oh_m91 gt -900.0),nm91_good)
    m91_good_upper = where((ohdust[mz[m91_good]].r23branch_m91 eq 'U'),nm91_good_upper)
    m91_good_lower = where((ohdust[mz[m91_good]].r23branch_m91 eq 'L'),nm91_good_lower)
    m91_good_ambig = where((ohdust[mz[m91_good]].r23branch_m91 eq 'A'),nm91_good_ambig)
    splog, 'M91 abundances: '+string(nm91_good,format='(I0)')+'/'+string(nmz,format='(I0)')+' ('+$
      string(round(100.0*nm91_good/nmz),format='(I0)')+'%).'
    splog, '   Upper: '+string(nm91_good_upper,format='(I0)')+'/'+string(nm91_good,format='(I0)')+' ('+$
      string(round(100.0*nm91_good_upper/nm91_good),format='(I0)')+'%).'
    splog, '   Lower: '+string(nm91_good_lower,format='(I0)')+'/'+string(nm91_good,format='(I0)')+' ('+$
      string(round(100.0*nm91_good_lower/nm91_good),format='(I0)')+'%).'
    splog, '   Ambig: '+string(nm91_good_ambig,format='(I0)')+'/'+string(nm91_good,format='(I0)')+' ('+$
      string(round(100.0*nm91_good_ambig/nm91_good),format='(I0)')+'%).'
    print

    kk04_good = where((ohdust[mz].zstrong_12oh_kk04 gt -900.0),nkk04_good)
    kk04_good_upper = where((ohdust[mz[kk04_good]].r23branch_kk04 eq 'U'),nkk04_good_upper)
    kk04_good_lower = where((ohdust[mz[kk04_good]].r23branch_kk04 eq 'L'),nkk04_good_lower)
    kk04_good_ambig = where((ohdust[mz[kk04_good]].r23branch_kk04 eq 'A'),nkk04_good_ambig)
    splog, 'KK04 abundances: '+string(nkk04_good,format='(I0)')+'/'+string(nmz,format='(I0)')+' ('+$
      string(round(100.0*nkk04_good/nmz),format='(I0)')+'%).'
    splog, '   Upper: '+string(nkk04_good_upper,format='(I0)')+'/'+string(nkk04_good,format='(I0)')+' ('+$
      string(round(100.0*nkk04_good_upper/nkk04_good),format='(I0)')+'%).'
    splog, '   Lower: '+string(nkk04_good_lower,format='(I0)')+'/'+string(nkk04_good,format='(I0)')+' ('+$
      string(round(100.0*nkk04_good_lower/nkk04_good),format='(I0)')+'%).'
    splog, '   Ambig: '+string(nkk04_good_ambig,format='(I0)')+'/'+string(nkk04_good,format='(I0)')+' ('+$
      string(round(100.0*nkk04_good_ambig/nkk04_good),format='(I0)')+'%).'
    print

; print some additional statistics

    goodmass = where((allsdssancillary[parent].kcorr_mass gt -900.0) and (allsdssancillary[parent].kauffmann_mass gt -900.0),ngoodmass)
    stats = im_stats(allsdssancillary[parent[goodmass]].kcorr_mass-allsdssancillary[parent[goodmass]].kauffmann_mass)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '+/-'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
    splog, 'Mass (k-correct minus Kauffmann) : '+xstr
;   plot, allsdssancillary[parent[goodmass]].kcorr_mass, allsdssancillary[parent[goodmass]].kauffmann_mass, ps=3, xr=[6,13], yr=[6,13]
;   oplot, !x.crange, !y.crange & oplot, !x.crange, !y.crange-stats.median_rej, line=2
    
    goodpt05 = where((ohnodust[mz].zstrong_12oh_pt05 gt -900.0) and (allsdssancillary[mz].tremonti_oh gt -900.0),ngoodpt05)
    stats = im_stats(ohnodust[mz[goodpt05]].zstrong_12oh_pt05-allsdssancillary[mz[goodpt05]].tremonti_oh)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '+/-'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
    splog, '12+log(O/H) (PT05 minus Tremonti): '+xstr
;   plot, ohnodust[mz[goodpt05]].zstrong_12oh_pt05, allsdssancillary[mz[goodpt05]].tremonti_oh, ps=3, xr=[7.0,9.6], yr=[7.0,9.6]
;   oplot, !x.crange, !y.crange & oplot, !x.crange, !y.crange-stats.median_rej, line=2
        
    goodm91 = where((ohnodust[mz].zstrong_12oh_m91 gt -900.0) and (allsdssancillary[mz].tremonti_oh gt -900.0),ngoodm91)
    stats = im_stats(ohnodust[mz[goodm91]].zstrong_12oh_m91-allsdssancillary[mz[goodm91]].tremonti_oh)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '+/-'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
    splog, '12+log(O/H) (M91 minus Tremonti) : '+xstr
;   plot, ohnodust[mz[goodm91]].zstrong_12oh_m91, allsdssancillary[mz[goodm91]].tremonti_oh, ps=3, xr=[7.0,9.6], yr=[7.0,9.6]
;   oplot, !x.crange, !y.crange & oplot, !x.crange, !y.crange-stats.median_rej, line=2
        
    goodkk04 = where((ohnodust[mz].zstrong_12oh_kk04 gt -900.0) and (allsdssancillary[mz].tremonti_oh gt -900.0),ngoodkk04)
    stats = im_stats(ohnodust[mz[goodkk04]].zstrong_12oh_kk04-allsdssancillary[mz[goodkk04]].tremonti_oh)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '+/-'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
    splog, '12+log(O/H) (KK04 minus Tremonti): '+xstr
;   plot, ohnodust[mz[goodkk04]].zstrong_12oh_kk04, allsdssancillary[mz[goodkk04]].tremonti_oh, ps=3, xr=[7.0,9.6], yr=[7.0,9.6]
;   oplot, !x.crange, !y.crange & oplot, !x.crange, !y.crange-stats.median_rej, line=2
        
    goodkk04 = where((ohnodust[mz].zstrong_12oh_kk04 gt -900.0) and (ohnodust[mz].zstrong_12oh_m91 gt -900.0),ngoodkk04)
    stats = im_stats(ohnodust[mz[goodkk04]].zstrong_12oh_kk04-ohnodust[mz[goodkk04]].zstrong_12oh_m91)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '+/-'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
    splog, '12+log(O/H) (KK04 minus M91)     : '+xstr
;   plot, ohnodust[mz[goodkk04]].zstrong_12oh_kk04, ohnodust[mz[goodkk04]].zstrong_12oh_m91, ps=3, xr=[7.0,9.6], yr=[7.0,9.6]
;   oplot, !x.crange, !y.crange & oplot, !x.crange, !y.crange-stats.median_rej, line=2
        
    goodkk04 = where((ohnodust[mz].zstrong_12oh_kk04 gt -900.0) and (ohnodust[mz].zstrong_12oh_pt05 gt -900.0),ngoodkk04)
    stats = im_stats(ohnodust[mz[goodkk04]].zstrong_12oh_kk04-ohnodust[mz[goodkk04]].zstrong_12oh_pt05)
    xstr = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '+/-'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'
    splog, '12+log(O/H) (KK04 minus PT05)    : '+xstr
;   plot, ohnodust[mz[goodkk04]].zstrong_12oh_kk04, ohnodust[mz[goodkk04]].zstrong_12oh_pt05, ps=3, xr=[7.0,9.6], yr=[7.0,9.6]
;   oplot, !x.crange, !y.crange & oplot, !x.crange, !y.crange-stats.median_rej, line=2
    
; write out; do some memory juggling
    
    if keyword_set(write) then begin

; ancillary       
       
       splog, 'Writing '+outpath+'sdss_abundances_hii_ancillary.fits.gz'
       mwrfits, allsdssancillary[mz[hii]], outpath+'sdss_abundances_hii_ancillary.fits', /create
       spawn, ['gzip -f '+outpath+'sdss_abundances_hii_ancillary.fits'], /sh

       splog, 'Writing '+outpath+'sdss_abundances_agn_ancillary.fits.gz'
       mwrfits, allsdssancillary[mz[agn]], outpath+'sdss_abundances_agn_ancillary.fits', /create
       spawn, ['gzip -f '+outpath+'sdss_abundances_agn_ancillary.fits'], /sh

       delvarx, allsdssancillary
       
; HII/AGN dust       

       splog, 'Appending additional tags.'
       allsdssdust = struct_addtags(temporary(allsdssdust),temporary(ohdust))

       splog, 'Writing '+outpath+'sdss_abundances_hii_speclinefit.fits.gz'
       mwrfits, allsdssdust[mz[hii]], outpath+'sdss_abundances_hii_speclinefit.fits', /create
       spawn, ['gzip -f '+outpath+'sdss_abundances_hii_speclinefit.fits'], /sh

       splog, 'Writing '+outpath+'sdss_abundances_agn_speclinefit.fits.gz'
       mwrfits, allsdssdust[mz[agn]], outpath+'sdss_abundances_agn_speclinefit.fits', /create
       spawn, ['gzip -f '+outpath+'sdss_abundances_agn_speclinefit.fits'], /sh

       delvarx, allsdssdust

; HII/AGN nodust       

       allsdssnodust = struct_addtags(temporary(allsdssnodust),temporary(ohnodust))

       splog, 'Writing '+outpath+'sdss_abundances_hii_speclinefit_nodust.fits.gz'
       mwrfits, allsdssnodust[mz[hii]], outpath+'sdss_abundances_hii_speclinefit_nodust.fits', /create
       spawn, ['gzip -f '+outpath+'sdss_abundances_hii_speclinefit_nodust.fits'], /sh

       splog, 'Writing '+outpath+'sdss_abundances_agn_speclinefit_nodust.fits.gz'
       mwrfits, allsdssnodust[mz[agn]], outpath+'sdss_abundances_agn_speclinefit_nodust.fits', /create
       spawn, ['gzip -f '+outpath+'sdss_abundances_agn_speclinefit_nodust.fits'], /sh

       delvarx, allsdssnodust

    endif

return
end
    
