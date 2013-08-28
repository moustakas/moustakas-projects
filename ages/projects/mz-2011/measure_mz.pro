function fit_mz, mass, oh, pivot=pivot, weight=weight, bin=bin, $
  minx=minx, maxx=maxx, minpts=minpts, verbose=verbose, $
  sample=sample, zlo=zlo, zhi=zhi, sdss=sdss, zoom=zoom, running=running
; fit the mass-metallicity relation
    
    number = n_elements(mass)
    if (n_elements(pivot) eq 0L) then pivot = 0.0
    if (n_elements(sample) eq 0L) then sample = ''
    if (n_elements(zlo) eq 0L) then zlo = 0.0
    if (n_elements(zhi) eq 0L) then zhi = 0.0
    
    if keyword_set(sdss) then $
      magrange = '$14.5<r_{AB}<17.6$' else $
      magrange = '$15<I<19.95$'

    running = im_medxbin(mass,oh,bin,minx=minx,maxx=maxx,$
      minpts=minpts,weight=weight,verbose=verbose)

    massaxis = findgen((max(running.binctr)-min(running.binctr))/$
      0.01+1)*0.01+min(running.binctr)
    coeff = poly_fit(running.binctr-pivot,running.medy,2)

    fit = {$
      sample:    sample, $
      zlo:       zlo,    $
      zhi:       zhi,    $
      magrange:  magrange,$
      coeff:     reform(coeff), $
;     coeff_err: 
      pivot:     pivot, $
      number:    number}

    if keyword_set(zoom) then begin
       xr=[9.7,11.5] & yr=[8.7,9.25]
    endif else begin
       xr=[7.8,11.6] & yr=[8.0,9.4]
    endelse
    hogg_scatterplot, mass, oh, ps=3, xsty=3, ysty=3, xr=xr, yr=yr, $
      xtitle=textoidl('log (M/M_{\odot})'), ytitle='12 + log (O/H)', $
      /internal
;      djs_plot, mass, oh, ps=3, xsty=3, ysty=3, xr=xr, yr=yr
    djs_oplot, running.binctr, running.medy, ps=7, color='red', sym=1.5, thick=6.0
    djs_oplot, massaxis, poly(massaxis-pivot,coeff), color='dark green', thick=8.0
    
    xchristy = [8.57,8.67,8.76,8.86,8.96,9.06,9.16,9.26,9.36,9.46,9.57,9.66,9.76,9.86,9.96,10.06,$
      10.16,10.26,10.36,10.46,10.56,10.66,10.76,10.86,10.95,11.05,11.15,11.25]
    ychristy = [8.44,8.48,8.57,8.61,8.63,8.66,8.68,8.71,8.74,8.78,8.82,8.84,8.87,8.90,8.94,$
      8.97,8.99,9.01,9.03,9.05,9.07,9.08,9.09,9.10,9.11,9.11,9.12,9.12]
    djs_oplot, xchristy, ychristy, psym=symcat(15), color='orange', symsize=1.0

    massaxis_t04 = findgen((11.3-8.5)/0.01+1)*0.01+8.5
    mzcoeff_t04 = [-1.492,1.847,-0.08026]
    mzcoeff_t04_shift = [-1.492,1.847,-0.08026]+[0.02,0.0,0.0]
    oh_t04 = poly(massaxis_t04,mzcoeff_t04)
    oh_t04_shift = poly(massaxis_t04,mzcoeff_t04_shift)
    oh_t04_salp = poly(massaxis_t04+im_convert_imf(/from_kroupa),mzcoeff_t04)
    djs_oplot, massaxis_t04, oh_t04, line=3, thick=6.0, color='blue'
;   djs_oplot, massaxis_t04, oh_t04_shift, line=3, thick=3.0, color='red'
;   djs_oplot, massaxis_t04, oh_t04_salp, line=3, thick=3.0, color='green'

return, fit
end

pro measure_mz, verbose=verbose

    sdsskcorr = read_sdss_mz_sample(/mzhii_kcorr)
    sdssohdust = read_sdss_mz_sample(/mzhii_log12oh)

    mpainfo = read_sdss_vagc_mpa(/postlss,sample='dr4',letter='all',poststr='0')
    mpaoh = read_sdss_vagc_mpa(/mpamassoh,sample='dr4',letter='all',poststr='0')

    spherematch, mpainfo.ra, mpainfo.dec, sdsskcorr.ra, $
      sdsskcorr.dec, 1.0/3600.0, m1, m2
    mpainfo1 = mpainfo[m1] & mpaoh1 = mpaoh[m1]
    sdsskcorr1 = sdsskcorr[m2] & sdssohdust1 = sdssohdust[m2]

    masspivot = 10.5
    minmass_mz_sdss = 9.4
    minmass_mz_sdss2 = 8.5
    maxmass_mz_sdss = 11.2
    minpts_mz_sdss = 100L
    binsize_mz_sdss = 0.1

    psfile = 'measure_mz.ps'
    im_plotfaves, /post & dfpsplot, psfile, /square, /color

; compare masses
    hogg_scatterplot, sdsskcorr1.mass+im_convert_imf(/from_chabrier), $
      mpaoh1.kauffmann_mass+im_convert_imf(/from_kroupa), $
      xr=[9.0,11.5], yr=[9.0,11.5], ytitle=textoidl('log (M/M_{\odot}) [Tremonti]'), $
      xtitle=textoidl('log (M/M_{\odot}) [Moustakas]'), /internal
    djs_oplot, !x.crange, !y.crange, line=0, thick=2.0
    massstats = im_stats((sdsskcorr1.mass+im_convert_imf(/from_chabrier))-$
      (mpaoh1.kauffmann_mass+im_convert_imf(/from_kroupa)),/verbose,sigrej=3.0)
; compare abundances
    ohgood = where((sdssohdust1.zstrong_ew_alpha_gr_12oh_kk04 gt -900.0) and $
      (mpaoh1.tremonti_oh gt -900.0))
    hogg_scatterplot, sdssohdust1[ohgood].zstrong_ew_alpha_gr_12oh_kk04, $
      mpaoh1[ohgood].tremonti_oh, $
      xr=[8.4,9.4], yr=[8.4,9.4], ytitle=textoidl('12 + log (O/H) [Tremonti]'), $
      xtitle=textoidl('12 + log (O/H) [Moustakas]'), /internal
    djs_oplot, !x.crange, !y.crange, line=0, thick=2.0
    ohstats = im_stats(mpaoh1[ohgood].tremonti_oh-sdssohdust1[ohgood].zstrong_ew_alpha_gr_12oh_kk04,/ver,sigrej=3.0)
; compare abundance residuals vs mass
    hogg_scatterplot, sdsskcorr1.mass+im_convert_imf(/from_chabrier), $
      mpaoh1.tremonti_oh-sdssohdust1.zstrong_ew_alpha_gr_12oh_kk04, $
      xr=[9.0,11.5], yr=[-0.25,0.25], ytitle=textoidl('\Delta [log (O/H)] [Tremonti - Moustakas]'), $
      xtitle=textoidl('log (M/M_{\odot}) [Moustakas]'), /internal
    djs_oplot, !x.crange, [0,0], line=0, thick=2.0

; Tremonti SDSS sample
    sdssindx = where((mpaoh.tremonti_oh gt -900.0) and (mpaoh.kauffmann_mass gt -900.0))
    sdssmass = mpaoh[sdssindx].kauffmann_mass; + im_convert_imf(/from_kroupa)
    sdssoh = mpaoh[sdssindx].tremonti_oh
    mz_sdss = fit_mz(sdssmass,sdssoh,pivot=masspivot,bin=binsize_mz_sdss,$
      minx=minmass_mz_sdss2,maxx=maxmass_mz_sdss,minpts=minpts_mz_sdss,$
      verbose=verbose,/sdss,running=t04)
; my SDSS sample    
    sdssindx = where((sdssohdust.zstrong_ew_alpha_gr_12oh_kk04 gt -900) and $
      (strtrim(sdssohdust.r23branch_ew_alpha_gr_kk04,2) eq 'U'))
    sdssmass = sdsskcorr[sdssindx].mass + im_convert_imf(/from_chabrier)
    sdssoh = sdssohdust[sdssindx].zstrong_ew_alpha_gr_12oh_kk04
; full range
    mz_sdss = fit_mz(sdssmass,sdssoh,pivot=masspivot,bin=binsize_mz_sdss,$
      minx=minmass_mz_sdss,maxx=maxmass_mz_sdss,minpts=minpts_mz_sdss,$
      verbose=verbose,/sdss,running=m08)
; zoom in
    mz_sdss = fit_mz(sdssmass,sdssoh,pivot=masspivot,bin=binsize_mz_sdss,$
      minx=minmass_mz_sdss,maxx=maxmass_mz_sdss,minpts=minpts_mz_sdss,$
      verbose=verbose,/sdss,/zoom)

; subtract the T04 and M08 median points
;   djs_plot, m08.binctr, m08.medy, ysty=3, ps=4, yr=[8.7,9.25]
;   djs_oplot, m08.binctr, interpol(t04.medy,t04.binctr,m08.binctr), ps=4, color='green'        
    djs_plot, m08.binctr, interpol(t04.medy,t04.binctr,m08.binctr)-$
      m08.medy+ohstats.median_rej*0.0, ysty=3, ps=-4, sym=2.5, $
      xtitle=textoidl('log (M/M_{\odot}) [Moustakas]'), $
      ytitle=textoidl('Median \Delta [log (O/H)] [Tremonti - Moustakas]')
    
    dfpsclose & im_plotfaves
    spawn, 'ps2pdf '+psfile+' ~/Desktop/'+repstr(psfile,'.ps','.pdf'), /sh


stop    
    
return
end
