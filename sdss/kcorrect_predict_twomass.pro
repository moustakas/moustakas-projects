pro zpt_hogg_scatterplot, mag, deltamag, band, input=input, _extra=extra
; make a scatterplot
    hogg_scatterplot, mag, deltamag, xsty=3, ysty=1, /internal, $
      xrange=[12.8,17.2], _extra=extra, $
      yrange=[-1,1], /outliers, $
      xtitle=textoidl(band+' (observed, AB mag)'), $
      ytitle=textoidl('\Delta'+band+' (observed minus synthesized, AB mag)')
    djs_oplot, !x.crange, [0,0], line=0, thick=5
    med = im_medxbin(mag,deltamag,0.2,minpts=10)
    djs_oplot, med.binctr, med.medy, line=0, thick=5, color='blue'
    djs_oplot, med.binctr, med.sigy84, line=5, thick=5, color='blue'
    djs_oplot, med.binctr, med.sigy16, line=5, thick=5, color='blue'
    im_legend, 'Input: '+input, /left, /top, box=0
    im_legend, '<\Delta'+band+'> = '+im_string_stats(deltamag,sigrej=3.0), $
      /left, /bottom, box=0

return
end

pro kcorrect_predict_twomass
; jm09may24nyu - predict and compare the 2MASS fluxes from just
;   fitting to the SDSS/ugriz (i.e., reproduce Fig. 10 in Blanton &
;   Roweis+07)

    common kcorrect_predict, phot1, tmass1, post1

    if (n_elements(post1) eq 0) then post1 = $
      mrdfits(vagc_name('post_catalog',sample='dr72',$
      letter='bsafe',post='25'),1)

    if (n_elements(phot1) eq 0) then phot1 = $
      hogg_mrdfits(vagc_name('object_sdss_imaging'),1,nrow=28800L)

    tfilters = 'twomass_'+['J','H','Ks']+'.par'
    if (n_elements(tmass1) eq 0) then tmass1 = $
      hogg_mrdfits(vagc_name('object_twomass'),1,nrow=28800L,$
      columns=['ra','decl','twomass_tag','k_m_ext','k_msig_ext',$
      'h_m_ext','h_msig_ext','j_m_ext','j_msig_ext'])

    post = post1
    tmass = tmass1[post.object_position]
    phot = phot1[post.object_position]
    keep = where((tmass.k_m_ext gt 0.0))
    keep = keep[0:10000]
    
    phot = phot[keep]
    tmass = tmass[keep]
    post = post[keep]

    twomass_to_maggies, tmass, tmaggies

; compute K-corrections using ugriz as input
    splog, 'K-correcting'
    sdss_to_maggies, inmaggies, inivarmaggies, calib=phot
    infilters = 'sdss_'+['u0','g0','r0','i0','z0']+'.par'

;   kcorr1 = im_kcorrect(post.z,inmaggies,inivarmaggies,$
;     infilters,band_shift=0.1,out_filterlist=tfilters,$
;     coeffs=coeffs1,synth_inmaggies=synth_inmaggies,$
;     synth_outmaggies=synth_outmaggies,/reset)

    t0 = systime(1)
    kcorr = sdss_kcorrect(post.z,calib=phot,band_shift=0.1,$
      coeffs=coeffs,rmaggies=rmaggies,mass=mass,$
      absmag=absmag,amivar=amivar,flux='petro')
    print, 'Time = ', systime(1)-t0

; synthesize 2MASS JHKs photometry and compare    
    k_load_vmatrix, vmatrix, lambda
    k_projection_table, rmatrix, vmatrix, lambda, zvals, $
      tfilters, zmin=0.01, zmax=0.25, nz=500L
    k_reconstruct_maggies, coeffs, post.z, synth_tmaggies, $
      rmatrix=rmatrix, zvals=zvals

    tband = ['J','H','Ks']
    magerr = [mean(tmass.j_msig_ext),mean(tmass.h_msig_ext),mean(tmass.k_msig_ext)]
    psfile = 'kcorrect_predict_twomass.ps'
    im_plotconfig, 0, pos, xmargin=[1.3,0.3], psfile=psfile, charsize=1.8
    for ii = 0, 2 do begin
       good = where((tmaggies[ii,*] gt 0.0) and (synth_tmaggies[ii,*] gt 0.0))
       mag = -2.5*alog10(tmaggies[ii,good])
       deltamag = -2.5*alog10(tmaggies[ii,good]/synth_tmaggies[ii,good])
       zpt_hogg_scatterplot, mag, deltamag, tband[ii], $
         input='ugriz', position=pos
       djs_oplot, !x.crange, +magerr[ii]*[1,1], line=5, thick=5.0
       djs_oplot, !x.crange, -magerr[ii]*[1,1], line=5, thick=5.0
    endfor
    im_plotconfig, psfile=psfile, /psclose, /gzip
    
stop    

return
end
