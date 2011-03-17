function t042dfunc, xy, p
    r23 = reform(xy[0,*])
    o32 = reform(xy[1,*])
    oh = poly(r23,p[0:3]) + poly(o32,p[4:5])
return, oh
end    

function fit_t042d, r23, o32, oh, oh_err=oh_err, $
  coeff_err=coeff_err, ohfit=ohfit, chi2=chi2

    nparams = 6
    parinfo = replicate({value: 1.0},nparams)

    xy = transpose([[r23],[o32]])
    coeff = mpfitfun('t042dfunc',xy,oh,oh_err,parinfo=parinfo,$
      perror=coeff_err,yfit=ohfit,bestnorm=bestnorm,quiet=0)
    chi2 = bestnorm/(n_elements(xy)-nparams)

return, coeff
end    

pro calibrate_t042d, ps=ps
; jm10aug12ucsd - calibrate a 2D version of the T04 calibration by
; fitting to R23 and O32 simultaneously 
;
; after much experimentation, it doesn't improve the
; metallicities that much!!

    ps = 1
    
    mzpath = ages_path(/projects)+'mz/'
    pspath = ages_path(/papers)+'mz/FIG_MZ/'
    if keyword_set(ps) then suffix = '.ps' else suffix = '.eps'

    sdss = read_mz_sample(/mzhii_ancillary,/sdss)
    ohnodust = read_mz_sample(/mzhii_log12oh,/nodust,/sdss)
    good = where((sdss.oh_median gt 0.0) and (ohnodust.r23 gt 0.0) and $
      (ohnodust.o32 gt 0.0),ngood)
    sdss = sdss[good]
    ohnodust = ohnodust[good]
    
    r23 = alog10(ohnodust.r23)
    o32 = alog10(ohnodust.o32)
    oh = sdss.oh_median
    oherr = (sdss.oh_p84-sdss.oh_p16)/2.0
    
; fit the model    
    r23axis = range(-0.1,0.9,500)
    t04_coeff = [9.185D,-0.313D,-0.264D,-0.321D]
;   coeff = fit_t042d(r23,o32,oh,oh_err=oh_err,ohfit=ohfit,chi2=chi2)
;   resid = oh-ohfit

;   coeff = poly_fit(r23,oh,3,measure_errors=oh_err,/double)
    im_poly_iter, r23, oh, 3, coeff=coeff, yerr=oh_err
;   coeff = robust_poly_fit(r23,oh,3,/double)
    resid = oh-poly(r23,coeff)
    rbin = im_medxbin(o32,resid,0.1,minpts=100,/verbose)

; make the plot    
    psfile = pspath+'calibrate_t042d'+suffix
    im_plotconfig, 6, pos, psfile=psfile, yspace=1.0, $
      xmargin=[1.3,0.4], width=6.8, height=[4.0,3.0]
; main plot
    mzsdss_hogg_scatterplot, r23, oh, position=pos[*,0], xsty=1, ysty=1, $
      xrange=[-0.3,1.1], yrange=[8,9.5], $
      xtitle=textoidl('log (R_{23})'), ytitle='12 + log (O/H)'
    djs_oplot, r23axis, poly(r23axis,t04_coeff), $
      line=5, color='red', thick=6
    djs_oplot, r23axis, poly(r23axis,coeff), $
      line=0, color='blue', thick=6
; residuals
    mzsdss_hogg_scatterplot, o32, resid, /noerase, position=pos[*,1], $
      xsty=1, ysty=1, xrange=[-1.1,0.5], yrange=0.6*[-1,1], $
      xtitle=textoidl('log (O_{32})'), ytitle='Residuals (dex)'
    oploterror, rbin.xbin, rbin.medy, rbin.sigy, psym=-symcat(6,thick=6), $
      symsize=3, errthick=6, color=fsc_color('navy',101), $
      errcolor=fsc_color('navy',101), thick=6
    im_plotconfig, /psclose
    
stop    

return
end
    
