pro sdss_bpt_archetypes, ll, sdss_rr
; jm09mar04nyu

    inpath = sdss_path(/mpa_dr7)
    outpath = sdss_path(/proj)+'archetype/'

;   ll = read_sdss_vagc_mpa(sample='dr7',poststr='32',/ispec)
    if (n_elements(ll) eq 0L) then $
      ll = mrdfits(inpath+'gal_line_dr7_v5_2.fit.gz',1)
    sdss = read_sdss_vagc_mpa(/ispec)

    spawn, 'grep \<variable lp_emline.cplex.out | '+$
      'grep value=\"1\" | cut -c 20-26', indx, /sh
    arch = parse_mpa_linefit(ll[indx])

    splog, 'Classifying'
    arch_cc = iclassification(arch,snrcut_class=0.0,ratios=arch_rr)
    ww = where((sdss.h_alpha[0]/sdss.h_alpha[1] gt 5.0) and $
      (sdss.h_beta[0]/sdss.h_beta[1] gt 5.0) and $
      (sdss.nii_6584[0]/sdss.nii_6584[1] gt 5.0) and $
      (sdss.oiii_5007[0]/sdss.oiii_5007[1] gt 5.0))
    if (n_elements(sdss_rr) eq 0L) then $
      sdss_cc = iclassification(sdss[ww],snrcut_class=5.0,ratios=sdss_rr)
    
    dfpsplot, 'sdss_bpt_archetypes.ps', /color, /square
    im_plotfaves, /post

    oiiihbrange = [-1.3,1.4]
    
; [NII]/Ha vs [OIII]/Hb    
    hogg_scatterplot, sdss_rr.nii_ha, sdss_rr.oiii_hb, /xsty, /ysty, $
      xrange=[-1.9,1.0], yrange=oiiihbrange, /outliers, outsymsize=0.05, $
      outcolor=djs_icolor('light blue'), $
      xtitle=textoidl('log ([N II] \lambda6584/H\alpha)'), $
      ytitle=textoidl('log ([O III] \lambda5007/H\beta)')
    djs_oplot, arch_rr.nii_ha, arch_rr.oiii_hb, psym=symcat(46,thick=5.0), $
      symsize=1.5, color='dark red'

    models = kewley_bpt_lines(/kauffmann,_extra=extra)
    oplot, models.x_nii, models.y_nii, line=0
    models = kewley_bpt_lines(_extra=extra)
    oplot, models.x_nii, models.y_nii, line=2

; [SII]/Ha vs [OIII]/Hb    
    hogg_scatterplot, alog10((sdss[ww].sii_6716[0]+sdss[ww].sii_6731[0])/sdss[ww].h_alpha[0]), $
      alog10(sdss[ww].oiii_5007[0]/sdss[ww].h_beta[0]), /xsty, /ysty, $
      xrange=[-1.5,0.5], yrange=oiiihbrange, /outliers, outsymsize=0.05, $
      outcolor=djs_icolor('light blue'), $
      xtitle=textoidl('log ([S II] \lambda\lambda6716,6731/H\alpha)'), $
      ytitle=textoidl('log ([O III] \lambda5007/H\beta)')
    djs_oplot, alog10((arch.sii_6716[0]+arch.sii_6731[0])/arch.h_alpha[0]), $
      alog10(arch.oiii_5007[0]/arch.h_beta[0]), psym=symcat(46,thick=5.0), $
      symsize=1.5, color='dark red'
    models = kewley_bpt_lines(_extra=extra)
    oplot, models.x_sii, models.y_sii, line=2
    
; [OI]/Ha vs [OIII]/Hb    
    hogg_scatterplot, alog10(sdss[ww].oi_6300[0]/sdss[ww].h_alpha[0]), $
      alog10(sdss[ww].oiii_5007[0]/sdss[ww].h_beta[0]), /xsty, /ysty, $
      xrange=[-2.5,0.0], yrange=oiiihbrange, /outliers, outsymsize=0.05, $
      outcolor=djs_icolor('light blue'), $
      xtitle=textoidl('log ([O I] \lambda6300/H\alpha)'), $
      ytitle=textoidl('log ([O III] \lambda5007/H\beta)')
    djs_oplot, alog10(arch.oi_6300[0]/arch.h_alpha[0]), $
      alog10(arch.oiii_5007[0]/arch.h_beta[0]), psym=symcat(46,thick=5.0), $
      symsize=1.5, color='dark red'
    models = kewley_bpt_lines(_extra=extra)
    oplot, models.x_oi, models.y_oi, line=2
    
    dfpsclose
    im_plotfaves

stop    
    
return
end
