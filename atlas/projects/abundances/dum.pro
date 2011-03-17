pro write_integrated_abundances_sample, atlasdust, atlasnodust, nfgsdust, nfgsnodust, write=write
; jm05aug26uofa - write integrated spectroscopic sample for the
;                 AGES/MZ paper 

    outpath = atlas_path(/projects)+'abundances/'
    snrcut = 3.0
    
    if (n_elements(atlasdust) eq 0L) then atlasdust = read_integrated(atlasnodust=atlasnodust)
    if (n_elements(nfgsdust) eq 0L) then nfgsdust = read_nfgs(nfgsnodust=nfgsnodust)

; which tags are in ATLASDUST but not in NFGSDUST

;   atlastags = tag_names(atlasdust[0]) & nfgstags = tag_names(nfgsdust[0])
;   niceprint, cmset_op(atlastags,'and',/not2,nfgstags)
;   atlastags = tag_names(atlasnodust[0]) & nfgstags = tag_names(nfgsnodust[0])
;   niceprint, cmset_op(atlastags,'and',/not2,nfgstags)

; append the two samples    
    
    intdust = struct_append(atlasdust,nfgsdust)
    intnodust = struct_append(atlasnodust,nfgsnodust)
    ngalaxy = n_elements(intdust)

; emission-line correction statistics

;   ugood = where((intdust.u_obs gt -900.0) and (intdust.u gt -900.0))
;   ustats = im_stats(intdust[ugood].u-intdust[ugood].u_obs,/verbose,/baremin)
;   bgood = where((intdust.b_obs gt -900.0) and (intdust.b gt -900.0))
;   bstats = im_stats(intdust[bgood].b-intdust[bgood].b_obs,/verbose,/baremin,/no_head)
;   vgood = where((intdust.v_obs gt -900.0) and (intdust.v gt -900.0))
;   vstats = im_stats(intdust[vgood].v-intdust[vgood].v_obs,/verbose,/baremin,/no_head)
;   rgood = where((intdust.r_obs gt -900.0) and (intdust.r gt -900.0))
;   rstats = im_stats(intdust[rgood].r-intdust[rgood].r_obs,/verbose,/baremin,/no_head)

; define more structure tags    
    
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

    mz = where((intdust.m_b_obs gt -900.0) and $
      (intdust.h_beta[0]/intdust.h_beta[1] gt snrcut) and $
      (intdust.h_alpha[0]/intdust.h_alpha[1] gt snrcut) and $
      (intdust.nii_6584[0]/intdust.nii_6584[1] gt snrcut) and $
      (intdust.oii_3727[0]/intdust.oii_3727[1] gt snrcut) and $
      (intdust.oiii_5007[0]/intdust.oiii_5007[1] gt snrcut),nmz)

    hii = where(strtrim(intdust[mz].bpt_pure_nii_class,2) eq 'HII',nhii)
    agn = where(strtrim(intdust[mz].bpt_pure_nii_class,2) eq 'AGN',nagn)
    splog, 'MZ sample: '+string(nmz,format='(I0)')+'/'+string(ngalaxy,format='(I0)')+' ('+$
      string(round(100.0*nmz/float(ngalaxy)),format='(I0)')+'%).'
    splog, 'MZ sample (star-forming): '+string(nhii,format='(I0)')+'/'+string(nmz,format='(I0)')+' ('+$
      string(round(100.0*nhii/float(nmz)),format='(I0)')+'%).'
    splog, 'MZ sample (type 2 AGN): '+string(nagn,format='(I0)')+'/'+string(nmz,format='(I0)')+' ('+$
      string(round(100.0*nagn/float(nmz)),format='(I0)')+'%).'
    print

; figure out the appropriate R23 branch; the EW abundances are the
; same for INTDUST and INTNODUST; for now, assign the
; unassigned objects to the upper branch, but check back later for an
; LZ method

; PT05/lower branch    
    
    lo_pt05_dust = where((intdust[mz].zstrong_niiha gt -900) and (intdust[mz].zstrong_niiha le -1.0) and $
                         (intdust[mz].zstrong_12oh_pt05_lower gt -900) and (intdust[mz].zstrong_12oh_pt05_upper gt -900) and $
                         (intdust[mz].zstrong_12oh_pt05_lower lt intdust[mz].zstrong_12oh_pt05_upper),nlo_pt05_dust)
    if (nlo_pt05_dust ne 0L) then begin
       ohdust[mz[lo_pt05_dust]].r23branch_pt05   = 'L'
       ohdust[mz[lo_pt05_dust]].zstrong_12oh_pt05          = intdust[mz[lo_pt05_dust]].zstrong_12oh_pt05_lower
       ohdust[mz[lo_pt05_dust]].zstrong_12oh_pt05_err      = intdust[mz[lo_pt05_dust]].zstrong_12oh_pt05_lower_err
    endif

    lo_pt05_nodust = where((intdust[mz].zstrong_niiha gt -900) and (intdust[mz].zstrong_niiha le -1.0) and $
                           (intnodust[mz].zstrong_12oh_pt05_lower gt -900) and (intnodust[mz].zstrong_12oh_pt05_upper gt -900) and $
                           (intnodust[mz].zstrong_12oh_pt05_lower lt intnodust[mz].zstrong_12oh_pt05_upper),nlo_pt05_nodust)
    if (nlo_pt05_nodust ne 0L) then begin
       ohnodust[mz[lo_pt05_nodust]].r23branch_pt05           = 'L'
       ohnodust[mz[lo_pt05_nodust]].zstrong_12oh_pt05        = intnodust[mz[lo_pt05_nodust]].zstrong_12oh_pt05_lower
       ohnodust[mz[lo_pt05_nodust]].zstrong_12oh_pt05_err    = intnodust[mz[lo_pt05_nodust]].zstrong_12oh_pt05_lower_err
    endif

    lo_pt05_ew = where((intdust[mz].zstrong_niiha gt -900) and (intdust[mz].zstrong_niiha le -1.0) and $
                       (intdust[mz].zstrong_ew_12oh_pt05_lower gt -900) and (intdust[mz].zstrong_ew_12oh_pt05_upper gt -900) and $
                       (intdust[mz].zstrong_ew_12oh_pt05_lower lt intdust[mz].zstrong_ew_12oh_pt05_upper),nlo_pt05_ew)
    if (nlo_pt05_ew ne 0L) then begin
       ohdust[mz[lo_pt05_ew]].r23branch_pt05_ew            = 'L'
       ohdust[mz[lo_pt05_ew]].zstrong_ew_12oh_pt05         = intdust[mz[lo_pt05_ew]].zstrong_ew_12oh_pt05_lower
       ohdust[mz[lo_pt05_ew]].zstrong_ew_12oh_pt05_err     = intdust[mz[lo_pt05_ew]].zstrong_ew_12oh_pt05_lower_err
       ohnodust[mz[lo_pt05_ew]].r23branch_pt05_ew          = 'L'
       ohnodust[mz[lo_pt05_ew]].zstrong_ew_12oh_pt05       = intdust[mz[lo_pt05_ew]].zstrong_ew_12oh_pt05_lower
       ohnodust[mz[lo_pt05_ew]].zstrong_ew_12oh_pt05_err   = intdust[mz[lo_pt05_ew]].zstrong_ew_12oh_pt05_lower_err
    endif

; PT05/upper branch    

    up_pt05_dust = where((intdust[mz].zstrong_niiha gt -900) and (intdust[mz].zstrong_niiha gt -1.0) and $
                         (intdust[mz].zstrong_12oh_pt05_upper gt -900) and (intdust[mz].zstrong_12oh_pt05_lower gt -900) and $
                         (intdust[mz].zstrong_12oh_pt05_upper gt intdust[mz].zstrong_12oh_pt05_lower),nup_pt05_dust)
    if (nup_pt05_dust ne 0L) then begin
       ohdust[mz[up_pt05_dust]].r23branch_pt05             = 'U'
       ohdust[mz[up_pt05_dust]].zstrong_12oh_pt05          = intdust[mz[up_pt05_dust]].zstrong_12oh_pt05_upper
       ohdust[mz[up_pt05_dust]].zstrong_12oh_pt05_err      = intdust[mz[up_pt05_dust]].zstrong_12oh_pt05_upper_err
    endif

    up_pt05_nodust = where((intdust[mz].zstrong_niiha gt -900) and (intdust[mz].zstrong_niiha gt -1.0) and $
                           (intnodust[mz].zstrong_12oh_pt05_upper gt -900) and (intnodust[mz].zstrong_12oh_pt05_lower gt -900) and $
                           (intnodust[mz].zstrong_12oh_pt05_upper gt intnodust[mz].zstrong_12oh_pt05_lower),nup_pt05_nodust)
    if (nup_pt05_nodust ne 0L) then begin
       ohnodust[mz[up_pt05_nodust]].r23branch_pt05           = 'U'
       ohnodust[mz[up_pt05_nodust]].zstrong_12oh_pt05        = intnodust[mz[up_pt05_nodust]].zstrong_12oh_pt05_upper
       ohnodust[mz[up_pt05_nodust]].zstrong_12oh_pt05_err    = intnodust[mz[up_pt05_nodust]].zstrong_12oh_pt05_upper_err
    endif

    up_pt05_ew = where((intdust[mz].zstrong_niiha gt -900) and (intdust[mz].zstrong_niiha gt -1.0) and $
                       (intdust[mz].zstrong_ew_12oh_pt05_lower gt -900) and (intdust[mz].zstrong_ew_12oh_pt05_upper gt -900) and $
                       (intdust[mz].zstrong_ew_12oh_pt05_upper gt intdust[mz].zstrong_ew_12oh_pt05_lower),nup_pt05_ew)
    if (nup_pt05_ew ne 0L) then begin
       ohdust[mz[up_pt05_ew]].r23branch_pt05_ew            = 'U'
       ohdust[mz[up_pt05_ew]].zstrong_ew_12oh_pt05         = intdust[mz[up_pt05_ew]].zstrong_ew_12oh_pt05_upper
       ohdust[mz[up_pt05_ew]].zstrong_ew_12oh_pt05_err     = intdust[mz[up_pt05_ew]].zstrong_ew_12oh_pt05_upper_err
       ohnodust[mz[up_pt05_ew]].r23branch_pt05_ew          = 'U'
       ohnodust[mz[up_pt05_ew]].zstrong_ew_12oh_pt05       = intdust[mz[up_pt05_ew]].zstrong_ew_12oh_pt05_upper
       ohnodust[mz[up_pt05_ew]].zstrong_ew_12oh_pt05_err   = intdust[mz[up_pt05_ew]].zstrong_ew_12oh_pt05_upper_err
    endif
    
; PT05/ambigiuos    

    ambig_pt05_dust = where((intdust[mz].zstrong_12oh_pt05_upper gt -900) and (intdust[mz].zstrong_12oh_pt05_lower gt -900) and $
                            (intdust[mz].zstrong_12oh_pt05_upper lt intdust[mz].zstrong_12oh_pt05_lower) and $
                            ((intdust[mz].zstrong_12oh_pt05_upper+intdust[mz].zstrong_12oh_pt05_upper_err) gt $
                             (intdust[mz].zstrong_12oh_pt05_lower-intdust[mz].zstrong_12oh_pt05_lower_err)),nambig_pt05_dust)

;   niceprint, (intdust[mz[ambig_pt05_dust]].zstrong_12oh_pt05_upper+intdust[mz[ambig_pt05_dust]].zstrong_12oh_pt05_upper_err), $
;     (intdust[mz[ambig_pt05_dust]].zstrong_12oh_pt05_lower-intdust[mz[ambig_pt05_dust]].zstrong_12oh_pt05_lower_err), $
;     intdust[mz[ambig_pt05_dust]].zstrong_12oh_pt05_upper_err, intdust[mz[ambig_pt05_dust]].zstrong_12oh_pt05_lower_err
;   plot, intdust[mz[ambig_pt05_dust]].zstrong_12oh_pt05_lower, intdust[mz[ambig_pt05_dust]].zstrong_12oh_pt05_upper, $
;     ps=4, xsty=3, ysty=3, xrange=[7.5,9.0], yrange=[7.5,9.0] & oplot, !x.crange, !y.crange
;   ploterror, intdust[mz[ambig_pt05_dust]].zstrong_12oh_pt05_lower, intdust[mz[ambig_pt05_dust]].zstrong_12oh_pt05_upper, $
;     intdust[mz[ambig_pt05_dust]].zstrong_12oh_pt05_lower_err, intdust[mz[ambig_pt05_dust]].zstrong_12oh_pt05_upper_err, $
;     ps=4, xsty=3, ysty=3, xrange=[7.5,9.0], yrange=[7.5,9.0], errstyle=1 & oplot, !x.crange, !y.crange

;   ambig_pt05_dust = where((intdust[mz].zstrong_12oh_pt05_upper gt -900) and (intdust[mz].zstrong_12oh_pt05_lower gt -900) and $
;                           ((intdust[mz].zstrong_12oh_pt05_upper lt intdust[mz].zstrong_12oh_pt05_lower) or $
;                            (intdust[mz].zstrong_12oh_pt05_lower gt intdust[mz].zstrong_12oh_pt05_upper)) and $
;                           (abs(intdust[mz].zstrong_12oh_pt05_upper-intdust[mz].zstrong_12oh_pt05_lower) lt deltaoh),nambig_pt05_dust)

    if (nambig_pt05_dust ne 0L) then begin
       ohdust[mz[ambig_pt05_dust]].r23branch_pt05 = 'A'
       for i = 0L, nambig_pt05_dust-1L do begin
          oh     = [intdust[mz[ambig_pt05_dust[i]]].zstrong_12oh_pt05_upper,$
                    intdust[mz[ambig_pt05_dust[i]]].zstrong_12oh_pt05_lower]
          oh_err = [intdust[mz[ambig_pt05_dust[i]]].zstrong_12oh_pt05_upper_err,$
                    intdust[mz[ambig_pt05_dust[i]]].zstrong_12oh_pt05_lower_err]
          ohdust[mz[ambig_pt05_dust[i]]].zstrong_12oh_pt05     = total(oh/oh_err^2)/total(1.0/oh_err^2)
          ohdust[mz[ambig_pt05_dust[i]]].zstrong_12oh_pt05_err = 1.0/sqrt(total(1.0/oh_err^2))
       endfor
    endif

;   ploterror, ohdust[mz[ambig_pt05_dust]].zstrong_12oh_pt05, intdust[mz[ambig_pt05_dust]].zstrong_12oh_pt05_lower, $
;     intdust[mz[ambig_pt05_dust]].zstrong_12oh_pt05_lower_err, ps=4, xsty=3, ysty=3, xrange=[7.3,8.6], yrange=[7.3,8.6], $
;     errstyle=1
;   oploterror, ohdust[mz[ambig_pt05_dust]].zstrong_12oh_pt05, intdust[mz[ambig_pt05_dust]].zstrong_12oh_pt05_upper, $
;     intdust[mz[ambig_pt05_dust]].zstrong_12oh_pt05_upper_err, ps=4, errstyle=1, color=djs_icolor('red'), errcolor=djs_icolor('red')
;   oplot, !x.crange, !y.crange
    
    ambig_pt05_nodust = where((intnodust[mz].zstrong_12oh_pt05_upper gt -900) and (intnodust[mz].zstrong_12oh_pt05_lower gt -900) and $
                              (intnodust[mz].zstrong_12oh_pt05_upper lt intnodust[mz].zstrong_12oh_pt05_lower) and $
                              ((intnodust[mz].zstrong_12oh_pt05_upper+intnodust[mz].zstrong_12oh_pt05_upper_err) gt $
                               (intnodust[mz].zstrong_12oh_pt05_lower-intnodust[mz].zstrong_12oh_pt05_lower_err)),nambig_pt05_nodust)

    if (nambig_pt05_nodust ne 0L) then begin
       ohnodust[mz[ambig_pt05_nodust]].r23branch_pt05 = 'A'
       for i = 0L, nambig_pt05_nodust-1L do begin
          oh     = [intnodust[mz[ambig_pt05_nodust[i]]].zstrong_12oh_pt05_upper,$
                    intnodust[mz[ambig_pt05_nodust[i]]].zstrong_12oh_pt05_lower]
          oh_err = [intnodust[mz[ambig_pt05_nodust[i]]].zstrong_12oh_pt05_upper_err,$
                    intnodust[mz[ambig_pt05_nodust[i]]].zstrong_12oh_pt05_lower_err]
          ohnodust[mz[ambig_pt05_nodust[i]]].zstrong_12oh_pt05     = total(oh/oh_err^2)/total(1.0/oh_err^2)
          ohnodust[mz[ambig_pt05_nodust[i]]].zstrong_12oh_pt05_err = 1.0/sqrt(total(1.0/oh_err^2))
       endfor
    endif

    ambig_pt05_ew = where((intdust[mz].zstrong_ew_12oh_pt05_upper gt -900) and (intdust[mz].zstrong_ew_12oh_pt05_lower gt -900) and $
                          (intdust[mz].zstrong_ew_12oh_pt05_upper lt intdust[mz].zstrong_ew_12oh_pt05_lower) and $
                          ((intdust[mz].zstrong_ew_12oh_pt05_upper+intdust[mz].zstrong_ew_12oh_pt05_upper_err) gt $
                           (intdust[mz].zstrong_ew_12oh_pt05_lower-intdust[mz].zstrong_ew_12oh_pt05_lower_err)),nambig_pt05_ew)

    if (nambig_pt05_ew ne 0L) then begin
       ohdust[mz[ambig_pt05_ew]].r23branch_pt05_ew   = 'A'
       ohnodust[mz[ambig_pt05_ew]].r23branch_pt05_ew = 'A'
       for i = 0L, nambig_pt05_ew-1L do begin
          oh     = [intdust[mz[ambig_pt05_ew[i]]].zstrong_ew_12oh_pt05_upper,$
                    intdust[mz[ambig_pt05_ew[i]]].zstrong_ew_12oh_pt05_lower]
          oh_err = [intdust[mz[ambig_pt05_ew[i]]].zstrong_ew_12oh_pt05_upper_err,$
                    intdust[mz[ambig_pt05_ew[i]]].zstrong_ew_12oh_pt05_lower_err]
          ohdust[mz[ambig_pt05_ew[i]]].zstrong_ew_12oh_pt05       = total(oh/oh_err^2)/total(1.0/oh_err^2)
          ohdust[mz[ambig_pt05_ew[i]]].zstrong_ew_12oh_pt05_err   = 1.0/sqrt(total(1.0/oh_err^2))
          ohnodust[mz[ambig_pt05_ew[i]]].zstrong_ew_12oh_pt05     = total(oh/oh_err^2)/total(1.0/oh_err^2)
          ohnodust[mz[ambig_pt05_ew[i]]].zstrong_ew_12oh_pt05_err = 1.0/sqrt(total(1.0/oh_err^2))
       endfor
    endif

; for now, assign everything else to the PT05 upper branch

    need_pt05_dust = where((intdust[mz].zstrong_12oh_pt05_upper gt -900.0) and (intdust[mz].zstrong_12oh_pt05_lower gt -900.0) and $
                           (intdust[mz].zstrong_12oh_pt05_upper gt intdust[mz].zstrong_12oh_pt05_lower) and $
                           (ohdust[mz].zstrong_12oh_pt05 lt -900.0),nneed_pt05_dust)
    if (nneed_pt05_dust ne 0L) then begin
       ohdust[mz[need_pt05_dust]].r23branch_pt05             = 'U'
       ohdust[mz[need_pt05_dust]].zstrong_12oh_pt05          = intdust[mz[need_pt05_dust]].zstrong_12oh_pt05_upper
       ohdust[mz[need_pt05_dust]].zstrong_12oh_pt05_err      = intdust[mz[need_pt05_dust]].zstrong_12oh_pt05_upper_err
    endif

    need_pt05_nodust = where((intnodust[mz].zstrong_12oh_pt05_upper gt -900.0) and (ohnodust[mz].zstrong_12oh_pt05 lt -900.0) and $
                             (intnodust[mz].zstrong_12oh_pt05_upper gt intnodust[mz].zstrong_12oh_pt05_lower) and $
                             (ohnodust[mz].zstrong_12oh_pt05 lt -900.0),nneed_pt05_nodust)
    if (nneed_pt05_nodust ne 0L) then begin
       ohnodust[mz[need_pt05_nodust]].r23branch_pt05             = 'U'
       ohnodust[mz[need_pt05_nodust]].zstrong_12oh_pt05          = intnodust[mz[need_pt05_nodust]].zstrong_12oh_pt05_upper
       ohnodust[mz[need_pt05_nodust]].zstrong_12oh_pt05_err      = intnodust[mz[need_pt05_nodust]].zstrong_12oh_pt05_upper_err
    endif

    need_pt05_ew_dust = where((intdust[mz].zstrong_ew_12oh_pt05_upper gt -900.0) and (intdust[mz].zstrong_ew_12oh_pt05_lower gt -900.0) and $
                              (intdust[mz].zstrong_ew_12oh_pt05_upper gt intdust[mz].zstrong_ew_12oh_pt05_lower) and $
                              (ohdust[mz].zstrong_ew_12oh_pt05 lt -900.0),nneed_pt05_ew_dust)
    if (nneed_pt05_ew_dust ne 0L) then begin
       ohdust[mz[need_pt05_ew_dust]].r23branch_pt05_ew        = 'U'
       ohdust[mz[need_pt05_ew_dust]].zstrong_ew_12oh_pt05     = intdust[mz[need_pt05_ew_dust]].zstrong_ew_12oh_pt05_upper
       ohdust[mz[need_pt05_ew_dust]].zstrong_ew_12oh_pt05_err = intdust[mz[need_pt05_ew_dust]].zstrong_ew_12oh_pt05_upper_err
    endif

    need_pt05_ew_nodust = where((intdust[mz].zstrong_ew_12oh_pt05_upper gt -900.0) and (intdust[mz].zstrong_ew_12oh_pt05_lower gt -900.0) and $
                                (intdust[mz].zstrong_ew_12oh_pt05_upper gt intdust[mz].zstrong_ew_12oh_pt05_lower) and $
                                (ohnodust[mz].zstrong_ew_12oh_pt05 lt -900.0),nneed_pt05_ew_nodust)
    if (nneed_pt05_ew_nodust ne 0L) then begin
       ohnodust[mz[need_pt05_ew_nodust]].r23branch_pt05_ew        = 'U'
       ohnodust[mz[need_pt05_ew_nodust]].zstrong_ew_12oh_pt05     = intdust[mz[need_pt05_ew_nodust]].zstrong_ew_12oh_pt05_upper
       ohnodust[mz[need_pt05_ew_nodust]].zstrong_ew_12oh_pt05_err = intdust[mz[need_pt05_ew_nodust]].zstrong_ew_12oh_pt05_upper_err
    endif

;   w = where((ohdust.zstrong_ew_12oh_pt05 gt -900) and (ohdust.r23branch_pt05_ew eq '?'),nw)
;   if (nw ne 0L) then stop
    
;   good = where((ohdust.zstrong_12oh_pt05 gt -900.0) and (ohdust.r23branch_pt05 eq '?'),ngood)
;   struct_print, struct_trimtags(ohdust[good],select='*PT05*')
;   if (ngood ne 0L) then stop

;   good = where(ohdust.zstrong_12oh_pt05 gt -900.0)
;   plot, intdust[good].zstrong_r23, ohdust[good].zstrong_12oh_pt05, ps=4, xsty=3, ysty=3
;   struct_print, struct_trimtags(ohdust[good],select='*PT05*')

; M91/lower branch    
    
    lo_m91_dust = where((intdust[mz].zstrong_niiha gt -900) and (intdust[mz].zstrong_niiha le -1.0) and $
                         (intdust[mz].zstrong_12oh_m91_lower gt -900) and (intdust[mz].zstrong_12oh_m91_upper gt -900) and $
                         (intdust[mz].zstrong_12oh_m91_lower lt intdust[mz].zstrong_12oh_m91_upper),nlo_m91_dust)
    if (nlo_m91_dust ne 0L) then begin
       ohdust[mz[lo_m91_dust]].r23branch_m91   = 'L'
       ohdust[mz[lo_m91_dust]].zstrong_12oh_m91          = intdust[mz[lo_m91_dust]].zstrong_12oh_m91_lower
       ohdust[mz[lo_m91_dust]].zstrong_12oh_m91_err      = intdust[mz[lo_m91_dust]].zstrong_12oh_m91_lower_err
    endif

    lo_m91_nodust = where((intdust[mz].zstrong_niiha gt -900) and (intdust[mz].zstrong_niiha le -1.0) and $
                           (intnodust[mz].zstrong_12oh_m91_lower gt -900) and (intnodust[mz].zstrong_12oh_m91_upper gt -900) and $
                           (intnodust[mz].zstrong_12oh_m91_lower lt intnodust[mz].zstrong_12oh_m91_upper),nlo_m91_nodust)
    if (nlo_m91_nodust ne 0L) then begin
       ohnodust[mz[lo_m91_nodust]].r23branch_m91           = 'L'
       ohnodust[mz[lo_m91_nodust]].zstrong_12oh_m91        = intnodust[mz[lo_m91_nodust]].zstrong_12oh_m91_lower
       ohnodust[mz[lo_m91_nodust]].zstrong_12oh_m91_err    = intnodust[mz[lo_m91_nodust]].zstrong_12oh_m91_lower_err
    endif

    lo_m91_ew = where((intdust[mz].zstrong_niiha gt -900) and (intdust[mz].zstrong_niiha le -1.0) and $
                       (intdust[mz].zstrong_ew_12oh_m91_lower gt -900) and (intdust[mz].zstrong_ew_12oh_m91_upper gt -900) and $
                       (intdust[mz].zstrong_ew_12oh_m91_lower lt intdust[mz].zstrong_ew_12oh_m91_upper),nlo_m91_ew)
    if (nlo_m91_ew ne 0L) then begin
       ohdust[mz[lo_m91_ew]].r23branch_m91_ew            = 'L'
       ohdust[mz[lo_m91_ew]].zstrong_ew_12oh_m91         = intdust[mz[lo_m91_ew]].zstrong_ew_12oh_m91_lower
       ohdust[mz[lo_m91_ew]].zstrong_ew_12oh_m91_err     = intdust[mz[lo_m91_ew]].zstrong_ew_12oh_m91_lower_err
       ohnodust[mz[lo_m91_ew]].r23branch_m91_ew          = 'L'
       ohnodust[mz[lo_m91_ew]].zstrong_ew_12oh_m91       = intdust[mz[lo_m91_ew]].zstrong_ew_12oh_m91_lower
       ohnodust[mz[lo_m91_ew]].zstrong_ew_12oh_m91_err   = intdust[mz[lo_m91_ew]].zstrong_ew_12oh_m91_lower_err
    endif

; M91/upper branch    

    up_m91_dust = where((intdust[mz].zstrong_niiha gt -900) and (intdust[mz].zstrong_niiha gt -1.0) and $
                         (intdust[mz].zstrong_12oh_m91_upper gt -900) and (intdust[mz].zstrong_12oh_m91_lower gt -900) and $
                         (intdust[mz].zstrong_12oh_m91_upper gt intdust[mz].zstrong_12oh_m91_lower),nup_m91_dust)
    if (nup_m91_dust ne 0L) then begin
       ohdust[mz[up_m91_dust]].r23branch_m91             = 'U'
       ohdust[mz[up_m91_dust]].zstrong_12oh_m91          = intdust[mz[up_m91_dust]].zstrong_12oh_m91_upper
       ohdust[mz[up_m91_dust]].zstrong_12oh_m91_err      = intdust[mz[up_m91_dust]].zstrong_12oh_m91_upper_err
    endif

    up_m91_nodust = where((intdust[mz].zstrong_niiha gt -900) and (intdust[mz].zstrong_niiha gt -1.0) and $
                           (intnodust[mz].zstrong_12oh_m91_upper gt -900) and (intnodust[mz].zstrong_12oh_m91_lower gt -900) and $
                           (intnodust[mz].zstrong_12oh_m91_upper gt intnodust[mz].zstrong_12oh_m91_lower),nup_m91_nodust)
    if (nup_m91_nodust ne 0L) then begin
       ohnodust[mz[up_m91_nodust]].r23branch_m91           = 'U'
       ohnodust[mz[up_m91_nodust]].zstrong_12oh_m91        = intnodust[mz[up_m91_nodust]].zstrong_12oh_m91_upper
       ohnodust[mz[up_m91_nodust]].zstrong_12oh_m91_err    = intnodust[mz[up_m91_nodust]].zstrong_12oh_m91_upper_err
    endif

    up_m91_ew = where((intdust[mz].zstrong_niiha gt -900) and (intdust[mz].zstrong_niiha gt -1.0) and $
                       (intdust[mz].zstrong_ew_12oh_m91_lower gt -900) and (intdust[mz].zstrong_ew_12oh_m91_upper gt -900) and $
                       (intdust[mz].zstrong_ew_12oh_m91_upper gt intdust[mz].zstrong_ew_12oh_m91_lower),nup_m91_ew)
    if (nup_m91_ew ne 0L) then begin
       ohdust[mz[up_m91_ew]].r23branch_m91_ew            = 'U'
       ohdust[mz[up_m91_ew]].zstrong_ew_12oh_m91         = intdust[mz[up_m91_ew]].zstrong_ew_12oh_m91_upper
       ohdust[mz[up_m91_ew]].zstrong_ew_12oh_m91_err     = intdust[mz[up_m91_ew]].zstrong_ew_12oh_m91_upper_err
       ohnodust[mz[up_m91_ew]].r23branch_m91_ew          = 'U'
       ohnodust[mz[up_m91_ew]].zstrong_ew_12oh_m91       = intdust[mz[up_m91_ew]].zstrong_ew_12oh_m91_upper
       ohnodust[mz[up_m91_ew]].zstrong_ew_12oh_m91_err   = intdust[mz[up_m91_ew]].zstrong_ew_12oh_m91_upper_err
    endif
    
; M91/ambigiuos    

    ambig_m91_dust = where((intdust[mz].zstrong_12oh_m91_upper gt -900) and (intdust[mz].zstrong_12oh_m91_lower gt -900) and $
                            (intdust[mz].zstrong_12oh_m91_upper lt intdust[mz].zstrong_12oh_m91_lower) and $
                            ((intdust[mz].zstrong_12oh_m91_upper+intdust[mz].zstrong_12oh_m91_upper_err) gt $
                             (intdust[mz].zstrong_12oh_m91_lower-intdust[mz].zstrong_12oh_m91_lower_err)),nambig_m91_dust)

    if (nambig_m91_dust ne 0L) then begin
       ohdust[mz[ambig_m91_dust]].r23branch_m91 = 'A'
       for i = 0L, nambig_m91_dust-1L do begin
          oh     = [intdust[mz[ambig_m91_dust[i]]].zstrong_12oh_m91_upper,$
                    intdust[mz[ambig_m91_dust[i]]].zstrong_12oh_m91_lower]
          oh_err = [intdust[mz[ambig_m91_dust[i]]].zstrong_12oh_m91_upper_err,$
                    intdust[mz[ambig_m91_dust[i]]].zstrong_12oh_m91_lower_err]
          ohdust[mz[ambig_m91_dust[i]]].zstrong_12oh_m91     = total(oh/oh_err^2)/total(1.0/oh_err^2)
          ohdust[mz[ambig_m91_dust[i]]].zstrong_12oh_m91_err = 1.0/sqrt(total(1.0/oh_err^2))
       endfor
    endif

    ambig_m91_nodust = where((intnodust[mz].zstrong_12oh_m91_upper gt -900) and (intnodust[mz].zstrong_12oh_m91_lower gt -900) and $
                             (intnodust[mz].zstrong_12oh_m91_upper lt intnodust[mz].zstrong_12oh_m91_lower) and $
                             ((intnodust[mz].zstrong_12oh_m91_upper+intnodust[mz].zstrong_12oh_m91_upper_err) gt $
                              (intnodust[mz].zstrong_12oh_m91_lower-intnodust[mz].zstrong_12oh_m91_lower_err)),nambig_m91_nodust)

    if (nambig_m91_nodust ne 0L) then begin
       ohnodust[mz[ambig_m91_nodust]].r23branch_m91 = 'A'
       for i = 0L, nambig_m91_nodust-1L do begin
          oh     = [intnodust[mz[ambig_m91_nodust[i]]].zstrong_12oh_m91_upper,$
                    intnodust[mz[ambig_m91_nodust[i]]].zstrong_12oh_m91_lower]
          oh_err = [intnodust[mz[ambig_m91_nodust[i]]].zstrong_12oh_m91_upper_err,$
                    intnodust[mz[ambig_m91_nodust[i]]].zstrong_12oh_m91_lower_err]
          ohnodust[mz[ambig_m91_nodust[i]]].zstrong_12oh_m91     = total(oh/oh_err^2)/total(1.0/oh_err^2)
          ohnodust[mz[ambig_m91_nodust[i]]].zstrong_12oh_m91_err = 1.0/sqrt(total(1.0/oh_err^2))
       endfor
    endif

    ambig_m91_ew = where((intdust[mz].zstrong_ew_12oh_m91_upper gt -900) and (intdust[mz].zstrong_ew_12oh_m91_lower gt -900) and $
                          (intdust[mz].zstrong_ew_12oh_m91_upper lt intdust[mz].zstrong_ew_12oh_m91_lower) and $
                          ((intdust[mz].zstrong_ew_12oh_m91_upper+intdust[mz].zstrong_ew_12oh_m91_upper_err) gt $
                           (intdust[mz].zstrong_ew_12oh_m91_lower-intdust[mz].zstrong_ew_12oh_m91_lower_err)),nambig_m91_ew)

    if (nambig_m91_ew ne 0L) then begin
       ohdust[mz[ambig_m91_ew]].r23branch_m91_ew   = 'A'
       ohnodust[mz[ambig_m91_ew]].r23branch_m91_ew = 'A'
       for i = 0L, nambig_m91_ew-1L do begin
          oh     = [intdust[mz[ambig_m91_ew[i]]].zstrong_ew_12oh_m91_upper,$
                    intdust[mz[ambig_m91_ew[i]]].zstrong_ew_12oh_m91_lower]
          oh_err = [intdust[mz[ambig_m91_ew[i]]].zstrong_ew_12oh_m91_upper_err,$
                    intdust[mz[ambig_m91_ew[i]]].zstrong_ew_12oh_m91_lower_err]
          ohdust[mz[ambig_m91_ew[i]]].zstrong_ew_12oh_m91       = total(oh/oh_err^2)/total(1.0/oh_err^2)
          ohdust[mz[ambig_m91_ew[i]]].zstrong_ew_12oh_m91_err   = 1.0/sqrt(total(1.0/oh_err^2))
          ohnodust[mz[ambig_m91_ew[i]]].zstrong_ew_12oh_m91     = total(oh/oh_err^2)/total(1.0/oh_err^2)
          ohnodust[mz[ambig_m91_ew[i]]].zstrong_ew_12oh_m91_err = 1.0/sqrt(total(1.0/oh_err^2))
       endfor
    endif

; for now, assign everything else to the M91 upper branch

    need_m91_dust = where((intdust[mz].zstrong_12oh_m91_upper gt -900.0) and (intdust[mz].zstrong_12oh_m91_lower gt -900.0) and $
                           (intdust[mz].zstrong_12oh_m91_upper gt intdust[mz].zstrong_12oh_m91_lower) and $
                           (ohdust[mz].zstrong_12oh_m91 lt -900.0),nneed_m91_dust)
    if (nneed_m91_dust ne 0L) then begin
       ohdust[mz[need_m91_dust]].r23branch_m91             = 'U'
       ohdust[mz[need_m91_dust]].zstrong_12oh_m91          = intdust[mz[need_m91_dust]].zstrong_12oh_m91_upper
       ohdust[mz[need_m91_dust]].zstrong_12oh_m91_err      = intdust[mz[need_m91_dust]].zstrong_12oh_m91_upper_err
    endif

    need_m91_nodust = where((intnodust[mz].zstrong_12oh_m91_upper gt -900.0) and (ohnodust[mz].zstrong_12oh_m91 lt -900.0) and $
                             (intnodust[mz].zstrong_12oh_m91_upper gt intnodust[mz].zstrong_12oh_m91_lower) and $
                             (ohnodust[mz].zstrong_12oh_m91 lt -900.0),nneed_m91_nodust)
    if (nneed_m91_nodust ne 0L) then begin
       ohnodust[mz[need_m91_nodust]].r23branch_m91             = 'U'
       ohnodust[mz[need_m91_nodust]].zstrong_12oh_m91          = intnodust[mz[need_m91_nodust]].zstrong_12oh_m91_upper
       ohnodust[mz[need_m91_nodust]].zstrong_12oh_m91_err      = intnodust[mz[need_m91_nodust]].zstrong_12oh_m91_upper_err
    endif

    need_m91_ew_dust = where((intdust[mz].zstrong_ew_12oh_m91_upper gt -900.0) and (intdust[mz].zstrong_ew_12oh_m91_lower gt -900.0) and $
                              (intdust[mz].zstrong_ew_12oh_m91_upper gt intdust[mz].zstrong_ew_12oh_m91_lower) and $
                              (ohdust[mz].zstrong_ew_12oh_m91 lt -900.0),nneed_m91_ew_dust)
    if (nneed_m91_ew_dust ne 0L) then begin
       ohdust[mz[need_m91_ew_dust]].r23branch_m91_ew        = 'U'
       ohdust[mz[need_m91_ew_dust]].zstrong_ew_12oh_m91     = intdust[mz[need_m91_ew_dust]].zstrong_ew_12oh_m91_upper
       ohdust[mz[need_m91_ew_dust]].zstrong_ew_12oh_m91_err = intdust[mz[need_m91_ew_dust]].zstrong_ew_12oh_m91_upper_err
    endif

    need_m91_ew_nodust = where((intdust[mz].zstrong_ew_12oh_m91_upper gt -900.0) and (intdust[mz].zstrong_ew_12oh_m91_lower gt -900.0) and $
                                (intdust[mz].zstrong_ew_12oh_m91_upper gt intdust[mz].zstrong_ew_12oh_m91_lower) and $
                                (ohnodust[mz].zstrong_ew_12oh_m91 lt -900.0),nneed_m91_ew_nodust)
    if (nneed_m91_ew_nodust ne 0L) then begin
       ohnodust[mz[need_m91_ew_nodust]].r23branch_m91_ew        = 'U'
       ohnodust[mz[need_m91_ew_nodust]].zstrong_ew_12oh_m91     = intdust[mz[need_m91_ew_nodust]].zstrong_ew_12oh_m91_upper
       ohnodust[mz[need_m91_ew_nodust]].zstrong_ew_12oh_m91_err = intdust[mz[need_m91_ew_nodust]].zstrong_ew_12oh_m91_upper_err
    endif

;   w = where((ohdust.zstrong_ew_12oh_m91 gt -900) and (ohdust.r23branch_m91_ew eq '?'),nw)
;   if (nw ne 0L) then stop
    
;   good = where((ohdust.zstrong_12oh_m91 gt -900.0) and (ohdust.r23branch_m91 eq '?'),ngood)
;   struct_print, struct_trimtags(ohdust[good],select='*M91*')
;   if (ngood ne 0L) then stop

;   good = where(ohdust.zstrong_12oh_m91 gt -900.0)
;   plot, intdust[good].zstrong_r23, ohdust[good].zstrong_12oh_m91, ps=4, xsty=3, ysty=3
;   struct_print, struct_trimtags(ohdust[good],select='*M91*')

; KK04/lower branch    
    
    lo_kk04_dust = where((intdust[mz].zstrong_niiha gt -900) and (intdust[mz].zstrong_niiha le -1.0) and $
                         (intdust[mz].zstrong_12oh_kk04_r23_lower gt -900) and (intdust[mz].zstrong_12oh_kk04_r23_upper gt -900) and $
                         (intdust[mz].zstrong_12oh_kk04_r23_lower lt intdust[mz].zstrong_12oh_kk04_r23_upper),nlo_kk04_dust)
    if (nlo_kk04_dust ne 0L) then begin
       ohdust[mz[lo_kk04_dust]].r23branch_kk04             = 'L'
       ohdust[mz[lo_kk04_dust]].zstrong_12oh_kk04          = intdust[mz[lo_kk04_dust]].zstrong_12oh_kk04_r23_lower
       ohdust[mz[lo_kk04_dust]].zstrong_12oh_kk04_err      = intdust[mz[lo_kk04_dust]].zstrong_12oh_kk04_r23_lower_err
    endif

    lo_kk04_nodust = where((intdust[mz].zstrong_niiha gt -900) and (intdust[mz].zstrong_niiha le -1.0) and $
                           (intnodust[mz].zstrong_12oh_kk04_r23_lower gt -900) and (intnodust[mz].zstrong_12oh_kk04_r23_upper gt -900) and $
                           (intnodust[mz].zstrong_12oh_kk04_r23_lower lt intnodust[mz].zstrong_12oh_kk04_r23_upper),nlo_kk04_nodust)
    if (nlo_kk04_nodust ne 0L) then begin
       ohnodust[mz[lo_kk04_nodust]].r23branch_kk04           = 'L'
       ohnodust[mz[lo_kk04_nodust]].zstrong_12oh_kk04        = intnodust[mz[lo_kk04_nodust]].zstrong_12oh_kk04_r23_lower
       ohnodust[mz[lo_kk04_nodust]].zstrong_12oh_kk04_err    = intnodust[mz[lo_kk04_nodust]].zstrong_12oh_kk04_r23_lower_err
    endif

    lo_kk04_ew = where((intdust[mz].zstrong_niiha gt -900) and (intdust[mz].zstrong_niiha le -1.0) and $
                       (intdust[mz].zstrong_ew_12oh_kk04_r23_lower gt -900) and (intdust[mz].zstrong_ew_12oh_kk04_r23_upper gt -900) and $
                       (intdust[mz].zstrong_ew_12oh_kk04_r23_lower lt intdust[mz].zstrong_ew_12oh_kk04_r23_upper),nlo_kk04_ew)
    if (nlo_kk04_ew ne 0L) then begin
       ohdust[mz[lo_kk04_ew]].r23branch_kk04_ew            = 'L'
       ohdust[mz[lo_kk04_ew]].zstrong_ew_12oh_kk04         = intdust[mz[lo_kk04_ew]].zstrong_ew_12oh_kk04_r23_lower
       ohdust[mz[lo_kk04_ew]].zstrong_ew_12oh_kk04_err     = intdust[mz[lo_kk04_ew]].zstrong_ew_12oh_kk04_r23_lower_err
       ohnodust[mz[lo_kk04_ew]].r23branch_kk04_ew          = 'L'
       ohnodust[mz[lo_kk04_ew]].zstrong_ew_12oh_kk04       = intdust[mz[lo_kk04_ew]].zstrong_ew_12oh_kk04_r23_lower
       ohnodust[mz[lo_kk04_ew]].zstrong_ew_12oh_kk04_err   = intdust[mz[lo_kk04_ew]].zstrong_ew_12oh_kk04_r23_lower_err
    endif
    
; KK04/upper branch    

    up_kk04_dust = where((intdust[mz].zstrong_niiha gt -900) and (intdust[mz].zstrong_niiha gt -1.0) and $
                         (intdust[mz].zstrong_12oh_kk04_r23_upper gt -900) and (intdust[mz].zstrong_12oh_kk04_r23_lower gt -900) and $
                         (intdust[mz].zstrong_12oh_kk04_r23_upper gt intdust[mz].zstrong_12oh_kk04_r23_lower),nup_kk04_dust)
    if (nup_kk04_dust ne 0L) then begin
       ohdust[mz[up_kk04_dust]].r23branch_kk04        = 'U'
       ohdust[mz[up_kk04_dust]].zstrong_12oh_kk04     = intdust[mz[up_kk04_dust]].zstrong_12oh_kk04_r23_upper
       ohdust[mz[up_kk04_dust]].zstrong_12oh_kk04_err = intdust[mz[up_kk04_dust]].zstrong_12oh_kk04_r23_upper_err
    endif

    up_kk04_nodust = where((intdust[mz].zstrong_niiha gt -900) and (intdust[mz].zstrong_niiha gt -1.0) and $
                           (intnodust[mz].zstrong_12oh_kk04_r23_upper gt -900) and (intnodust[mz].zstrong_12oh_kk04_r23_lower gt -900) and $
                           (intnodust[mz].zstrong_12oh_kk04_r23_upper gt intnodust[mz].zstrong_12oh_kk04_r23_lower),nup_kk04_nodust)
    if (nup_kk04_nodust ne 0L) then begin
       ohnodust[mz[up_kk04_nodust]].r23branch_kk04        = 'U'
       ohnodust[mz[up_kk04_nodust]].zstrong_12oh_kk04     = intnodust[mz[up_kk04_nodust]].zstrong_12oh_kk04_r23_upper
       ohnodust[mz[up_kk04_nodust]].zstrong_12oh_kk04_err = intnodust[mz[up_kk04_nodust]].zstrong_12oh_kk04_r23_upper_err
    endif

    up_kk04_ew = where((intdust[mz].zstrong_niiha gt -900) and (intdust[mz].zstrong_niiha gt -1.0) and $
                       (intdust[mz].zstrong_ew_12oh_kk04_r23_lower gt -900) and (intdust[mz].zstrong_ew_12oh_kk04_r23_upper gt -900) and $
                       (intdust[mz].zstrong_ew_12oh_kk04_r23_upper gt intdust[mz].zstrong_ew_12oh_kk04_r23_lower),nup_kk04_ew)
    if (nup_kk04_ew ne 0L) then begin
       ohdust[mz[up_kk04_ew]].r23branch_kk04_ew          = 'U'
       ohdust[mz[up_kk04_ew]].zstrong_ew_12oh_kk04       = intdust[mz[up_kk04_ew]].zstrong_ew_12oh_kk04_r23_upper
       ohdust[mz[up_kk04_ew]].zstrong_ew_12oh_kk04_err   = intdust[mz[up_kk04_ew]].zstrong_ew_12oh_kk04_r23_upper_err
       ohnodust[mz[up_kk04_ew]].r23branch_kk04_ew        = 'U'
       ohnodust[mz[up_kk04_ew]].zstrong_ew_12oh_kk04     = intdust[mz[up_kk04_ew]].zstrong_ew_12oh_kk04_r23_upper
       ohnodust[mz[up_kk04_ew]].zstrong_ew_12oh_kk04_err = intdust[mz[up_kk04_ew]].zstrong_ew_12oh_kk04_r23_upper_err
    endif

; KK04/ambigiuos    

    ambig_kk04_dust = where((intdust[mz].zstrong_12oh_kk04_r23_upper gt -900) and (intdust[mz].zstrong_12oh_kk04_r23_lower gt -900) and $
                            (intdust[mz].zstrong_12oh_kk04_r23_upper lt intdust[mz].zstrong_12oh_kk04_r23_lower) and $
                            ((intdust[mz].zstrong_12oh_kk04_r23_upper+intdust[mz].zstrong_12oh_kk04_r23_upper_err) gt $
                             (intdust[mz].zstrong_12oh_kk04_r23_lower-intdust[mz].zstrong_12oh_kk04_r23_lower_err)),nambig_kk04_dust)

    if (nambig_kk04_dust ne 0L) then begin
       ohdust[mz[ambig_kk04_dust]].r23branch_kk04 = 'A'
       for i = 0L, nambig_kk04_dust-1L do begin
          oh     = [intdust[mz[ambig_kk04_dust[i]]].zstrong_12oh_kk04_r23_upper,$
                    intdust[mz[ambig_kk04_dust[i]]].zstrong_12oh_kk04_r23_lower]
          oh_err = [intdust[mz[ambig_kk04_dust[i]]].zstrong_12oh_kk04_r23_upper_err,$
                    intdust[mz[ambig_kk04_dust[i]]].zstrong_12oh_kk04_r23_lower_err]
          ohdust[mz[ambig_kk04_dust[i]]].zstrong_12oh_kk04     = total(oh/oh_err^2)/total(1.0/oh_err^2)
          ohdust[mz[ambig_kk04_dust[i]]].zstrong_12oh_kk04_err = 1.0/sqrt(total(1.0/oh_err^2))
       endfor
    endif

    ambig_kk04_nodust = where((intnodust[mz].zstrong_12oh_kk04_r23_upper gt -900) and (intnodust[mz].zstrong_12oh_kk04_r23_lower gt -900) and $
                            (intnodust[mz].zstrong_12oh_kk04_r23_upper lt intnodust[mz].zstrong_12oh_kk04_r23_lower) and $
                            ((intnodust[mz].zstrong_12oh_kk04_r23_upper+intnodust[mz].zstrong_12oh_kk04_r23_upper_err) gt $
                             (intnodust[mz].zstrong_12oh_kk04_r23_lower-intnodust[mz].zstrong_12oh_kk04_r23_lower_err)),nambig_kk04_nodust)

    if (nambig_kk04_nodust ne 0L) then begin
       ohnodust[mz[ambig_kk04_nodust]].r23branch_kk04 = 'A'
       for i = 0L, nambig_kk04_nodust-1L do begin
          oh     = [intnodust[mz[ambig_kk04_nodust[i]]].zstrong_12oh_kk04_r23_upper,$
                    intnodust[mz[ambig_kk04_nodust[i]]].zstrong_12oh_kk04_r23_lower]
          oh_err = [intnodust[mz[ambig_kk04_nodust[i]]].zstrong_12oh_kk04_r23_upper_err,$
                    intnodust[mz[ambig_kk04_nodust[i]]].zstrong_12oh_kk04_r23_lower_err]
          ohnodust[mz[ambig_kk04_nodust[i]]].zstrong_12oh_kk04     = total(oh/oh_err^2)/total(1.0/oh_err^2)
          ohnodust[mz[ambig_kk04_nodust[i]]].zstrong_12oh_kk04_err = 1.0/sqrt(total(1.0/oh_err^2))
       endfor
    endif

    ambig_kk04_ew = where((intdust[mz].zstrong_ew_12oh_kk04_r23_upper gt -900) and (intdust[mz].zstrong_ew_12oh_kk04_r23_lower gt -900) and $
                          (intdust[mz].zstrong_ew_12oh_kk04_r23_upper lt intdust[mz].zstrong_ew_12oh_kk04_r23_lower) and $
                          ((intdust[mz].zstrong_ew_12oh_kk04_r23_upper+intdust[mz].zstrong_ew_12oh_kk04_r23_upper_err) gt $
                           (intdust[mz].zstrong_ew_12oh_kk04_r23_lower-intdust[mz].zstrong_ew_12oh_kk04_r23_lower_err)),nambig_kk04_ew)

    if (nambig_kk04_ew ne 0L) then begin
       ohdust[mz[ambig_kk04_ew]].r23branch_kk04_ew   = 'A'
       ohnodust[mz[ambig_kk04_ew]].r23branch_kk04_ew = 'A'
       for i = 0L, nambig_kk04_ew-1L do begin
          oh     = [intdust[mz[ambig_kk04_ew[i]]].zstrong_ew_12oh_kk04_r23_upper,$
                    intdust[mz[ambig_kk04_ew[i]]].zstrong_ew_12oh_kk04_r23_lower]
          oh_err = [intdust[mz[ambig_kk04_ew[i]]].zstrong_ew_12oh_kk04_r23_upper_err,$
                    intdust[mz[ambig_kk04_ew[i]]].zstrong_ew_12oh_kk04_r23_lower_err]
          ohdust[mz[ambig_kk04_ew[i]]].zstrong_ew_12oh_kk04       = total(oh/oh_err^2)/total(1.0/oh_err^2)
          ohdust[mz[ambig_kk04_ew[i]]].zstrong_ew_12oh_kk04_err   = 1.0/sqrt(total(1.0/oh_err^2))
          ohnodust[mz[ambig_kk04_ew[i]]].zstrong_ew_12oh_kk04     = total(oh/oh_err^2)/total(1.0/oh_err^2)
          ohnodust[mz[ambig_kk04_ew[i]]].zstrong_ew_12oh_kk04_err = 1.0/sqrt(total(1.0/oh_err^2))
       endfor
    endif

; for now, assign everything else to the upper branch

    need_kk04_dust = where((intdust[mz].zstrong_12oh_kk04_r23_upper gt -900.0) and (intdust[mz].zstrong_12oh_kk04_r23_lower gt -900.0) and $
                           (intdust[mz].zstrong_12oh_kk04_r23_upper gt intdust[mz].zstrong_12oh_kk04_r23_lower) and $
                           (ohdust[mz].zstrong_12oh_kk04 lt -900.0),nneed_kk04_dust)
    if (nneed_kk04_dust ne 0L) then begin
       ohdust[mz[need_kk04_dust]].r23branch_kk04        = 'U'
       ohdust[mz[need_kk04_dust]].zstrong_12oh_kk04     = intdust[mz[need_kk04_dust]].zstrong_12oh_kk04_r23_upper
       ohdust[mz[need_kk04_dust]].zstrong_12oh_kk04_err = intdust[mz[need_kk04_dust]].zstrong_12oh_kk04_r23_upper_err
    endif

    need_kk04_nodust = where((intnodust[mz].zstrong_12oh_kk04_r23_upper gt -900.0) and (ohnodust[mz].zstrong_12oh_kk04 lt -900.0) and $
                             (intnodust[mz].zstrong_12oh_kk04_r23_upper gt intnodust[mz].zstrong_12oh_kk04_r23_lower) and $
                             (ohnodust[mz].zstrong_12oh_kk04 lt -900.0),nneed_kk04_nodust)
    if (nneed_kk04_nodust ne 0L) then begin
       ohnodust[mz[need_kk04_nodust]].r23branch_kk04        = 'U'
       ohnodust[mz[need_kk04_nodust]].zstrong_12oh_kk04     = intnodust[mz[need_kk04_nodust]].zstrong_12oh_kk04_r23_upper
       ohnodust[mz[need_kk04_nodust]].zstrong_12oh_kk04_err = intnodust[mz[need_kk04_nodust]].zstrong_12oh_kk04_r23_upper_err
    endif

    need_kk04_ew_dust = where((intdust[mz].zstrong_ew_12oh_kk04_r23_upper gt -900.0) and (intdust[mz].zstrong_ew_12oh_kk04_r23_lower gt -900.0) and $
                              (intdust[mz].zstrong_ew_12oh_kk04_r23_upper gt intdust[mz].zstrong_ew_12oh_kk04_r23_lower) and $
                              (ohdust[mz].zstrong_ew_12oh_kk04 lt -900.0),nneed_kk04_ew_dust)
    if (nneed_kk04_ew_dust ne 0L) then begin
       ohdust[mz[need_kk04_ew_dust]].r23branch_kk04_ew        = 'U'
       ohdust[mz[need_kk04_ew_dust]].zstrong_ew_12oh_kk04     = intdust[mz[need_kk04_ew_dust]].zstrong_ew_12oh_kk04_r23_upper
       ohdust[mz[need_kk04_ew_dust]].zstrong_ew_12oh_kk04_err = intdust[mz[need_kk04_ew_dust]].zstrong_ew_12oh_kk04_r23_upper_err
    endif

    need_kk04_ew_nodust = where((intdust[mz].zstrong_ew_12oh_kk04_r23_upper gt -900.0) and (intdust[mz].zstrong_ew_12oh_kk04_r23_lower gt -900.0) and $
                              (intdust[mz].zstrong_ew_12oh_kk04_r23_upper gt intdust[mz].zstrong_ew_12oh_kk04_r23_lower) and $
                              (ohnodust[mz].zstrong_ew_12oh_kk04 lt -900.0),nneed_kk04_ew_nodust)
    if (nneed_kk04_ew_nodust ne 0L) then begin
       ohnodust[mz[need_kk04_ew_nodust]].r23branch_kk04_ew        = 'U'
       ohnodust[mz[need_kk04_ew_nodust]].zstrong_ew_12oh_kk04     = intdust[mz[need_kk04_ew_nodust]].zstrong_ew_12oh_kk04_r23_upper
       ohnodust[mz[need_kk04_ew_nodust]].zstrong_ew_12oh_kk04_err = intdust[mz[need_kk04_ew_nodust]].zstrong_ew_12oh_kk04_r23_upper_err
    endif

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

    if keyword_set(write) then begin

       splog, 'Appending additional tags.'
       intdust = struct_addtags(temporary(intdust),struct_addtags(moretags,ohdust))
       intnodust = struct_addtags(temporary(intnodust),struct_addtags(moretags,ohnodust))

; --------------------
       
       splog, 'Writing '+outpath+'integrated_abundances_speclinefit.fits.gz'
       mwrfits, intdust[mz], outpath+'integrated_abundances_speclinefit.fits', /create
       spawn, ['gzip -f '+outpath+'integrated_abundances_speclinefit.fits'], /sh

       splog, 'Writing '+outpath+'integrated_abundances_speclinefit_nodust.fits.gz'
       mwrfits, intnodust[mz], outpath+'integrated_abundances_speclinefit_nodust.fits', /create
       spawn, ['gzip -f '+outpath+'integrated_abundances_speclinefit_nodust.fits'], /sh

; --------------------
       
       splog, 'Writing '+outpath+'integrated_abundances_hii_speclinefit.fits.gz'
       mwrfits, intdust[mz[hii]], outpath+'integrated_abundances_hii_speclinefit.fits', /create
       spawn, ['gzip -f '+outpath+'integrated_abundances_hii_speclinefit.fits'], /sh

       splog, 'Writing '+outpath+'integrated_abundances_hii_speclinefit_nodust.fits.gz'
       mwrfits, intnodust[mz[hii]], outpath+'integrated_abundances_hii_speclinefit_nodust.fits', /create
       spawn, ['gzip -f '+outpath+'integrated_abundances_hii_speclinefit_nodust.fits'], /sh

; --------------------

       splog, 'Writing '+outpath+'integrated_abundances_agn_speclinefit.fits.gz'
       mwrfits, intdust[mz[agn]], outpath+'integrated_abundances_agn_speclinefit.fits', /create
       spawn, ['gzip -f '+outpath+'integrated_abundances_agn_speclinefit.fits'], /sh

       splog, 'Writing '+outpath+'integrated_abundances_agn_speclinefit_nodust.fits.gz'
       mwrfits, intnodust[mz[agn]], outpath+'integrated_abundances_agn_speclinefit_nodust.fits', /create
       spawn, ['gzip -f '+outpath+'integrated_abundances_agn_speclinefit_nodust.fits'], /sh

    endif

return
end
    
