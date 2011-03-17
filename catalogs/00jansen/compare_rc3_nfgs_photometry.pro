pro compare_rc3_nfgs_photometry, nfgs, postscript=postscript
; jm06jan16uofa - compare Rolf's extrapolated (total) B-band
;                 magnitudes against the RC3 photometry

    if (n_elements(nfgs) eq 0L) then nfgs = read_00jansen()
    nfgs_obs = nfgs

    pspath = '/home/ioannis/catalogs/00jansen/'
    
    if keyword_set(postscript) then begin
       postthick = 5.0 
    endif else begin
       postthick = 2.0 
    endelse
    
    pagemaker, nx=1, ny=2, xspace=0, yspace=0.0, width=6.5, height=[5.0,2.5], $
      xmargin=[1.5,0.5], ymargin=[0.4,1.1], xpage=8.5, ypage=9.0, $
      position=pos, /normal

    bmagrange = [11.7,16]
    
; ---------------------------------------------------------------------------
; B_T - "corrected"
; ---------------------------------------------------------------------------

    psname = 'nfgs_rc3_b'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=9.0, encapsulated=encapsulated

    indx_bt = where(strmatch(nfgs.rc3_b_origin,'*B_T*',/fold))
    nfgs_bt = nfgs[indx_bt].b
    rc3_bt = nfgs[indx_bt].rc3_b
    resid_bt = rc3_bt - nfgs_bt

    nfgs_obs[indx_bt].b = nfgs[indx_bt].b + nfgs[indx_bt].ebv_mw*k_lambda(4344.0,/odonnell)
    nfgs_obs[indx_bt].rc3_b = nfgs[indx_bt].rc3_b + nfgs[indx_bt].ebv_mw*k_lambda(4344.0,/odonnell)

    indx_mb = where(strmatch(nfgs.rc3_b_origin,'*m_B*',/fold))
    nfgs_mb = nfgs[indx_mb].b
    rc3_mb = nfgs[indx_mb].rc3_b
    resid_mb = rc3_mb - nfgs_mb

    nfgs_obs[indx_mb].b = nfgs[indx_mb].b + nfgs[indx_mb].ebv_mw*k_lambda(4344.0,/odonnell)
    nfgs_obs[indx_mb].rc3_b = nfgs[indx_mb].rc3_b + nfgs[indx_mb].ebv_mw*k_lambda(4344.0,/odonnell)

    stats = im_stats(resid_bt,/verbose)
    xstr_bt = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'

    stats = im_stats(resid_mb,/verbose)
    xstr_mb = strtrim(string(stats.median_rej,format='(F12.2)'),2)+' ('+$
      strtrim(string(stats.mean_rej,format='(F12.2)'),2)+$
      '\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)+')'

    djs_plot, [0], [0], /nodata, ps=8, xsty=3, ysty=3, charsize=1.7, charthick=postthick, $
      xthick=postthick, ythick=postthick, xrange=bmagrange, yrange=bmagrange, $
      position=pos[*,0], ytitle='B_{RC3} [Extinction-Corrected, mag]', xtitle='', xtickname=replicate(' ',10)
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick

    plotsym, 8, 1, fill=0, color=djs_icolor('blue')
    djs_oplot, nfgs_mb, rc3_mb, ps=8
    plotsym, 0, 1, /fill, color=djs_icolor('red')
    djs_oplot, nfgs_bt, rc3_bt, ps=8

    resrange = (max(abs(resid_bt)))*[-1.1,1.1]
    
    djs_plot, [0], [0], /nodata, xsty=3, ysty=3, ps=8, xrange=bmagrange, $
      position=pos[*,1], /noerase, yrange=resrange, charsize=1.7, $
      charthick=postthick, xthick=postthick, ythick=postthick, $
      xtitle='B_{tot, NFGS} [Extinction-Corrected, mag]', ytitle='Residuals [mag]'
    djs_oplot, !x.crange, [0,0], line=0, thick=postthick

    plotsym, 8, 1, fill=0, color=djs_icolor('blue')
    djs_oplot, nfgs_mb, resid_mb, ps=8
    plotsym, 0, 1, /fill, color=djs_icolor('red')
    djs_oplot, nfgs_bt, resid_bt, ps=8

    legend, textoidl(['B_{T}: '+xstr_bt,'m_{B}: '+xstr_mb]), /right, /top, box=0, charsize=1.3, $
      charthick=postthick, textcolor=djs_icolor(['red','blue'])

; output ascii file

    indx = where((strmatch(nfgs.rc3_b_origin,'*B_T*',/fold)))
;   indx = where((strmatch(nfgs.rc3_b_origin,'*B_T*',/fold)) or $
;     (strmatch(nfgs.rc3_b_origin,'*m_B*',/fold)),nindx)
    if keyword_set(postscript) then openw, lun, 'nfgs_rc3_b_cor.dat', /get_lun
    struct_print, struct_trimtags(nfgs[indx],select=['NFGS_ID','NFGS_GALAXY',$
      'B','B_ERR','RC3_B','RC3_B_ERR','RC3_B_ORIGIN']), lun=lun
    if keyword_set(postscript) then free_lun, lun
    
    if keyword_set(postscript) then openw, lun, 'nfgs_rc3_b_obs.dat', /get_lun
    struct_print, struct_trimtags(nfgs_obs[indx],select=['NFGS_ID','NFGS_GALAXY',$
      'B','B_ERR','RC3_B','RC3_B_ERR','RC3_B_ORIGIN']), lun=lun
    if keyword_set(postscript) then free_lun, lun
    
    im_openclose, postscript=postscript, /close
    
stop    
    
return
end
    
