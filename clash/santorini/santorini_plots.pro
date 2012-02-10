pro santorini_plots, supergrid, models=models, isedfit=isedfit, $
  qaplot=qaplot, clobber=clobber
; jm11nov08ucsd - build plots for the paper

;   common com_post, model, ised, mstar, post, age, Z, tau, sfr0, b100, av, sfrage
    
    datapath = clash_path(/santorini)
    isedpath = datapath+'isedfit/'

    isedfit_sfhgrid_dir = datapath+'montegrids/'
    sfhgrid_paramfile = getenv('CLASH_DIR')+'/santorini/santorini_sfhgrid.par'

    filters = clash_filterlist(nice=nice_filters,/useirac,weff=weff)
    ndraw = isedfit_ndraw()

;; gather the photometry and restore the results
;    cat = read_santorini()
;    nobj = n_elements(cat)

;   if (n_elements(mstar) eq 0L) then begin
;      mstar = isedfit_reconstruct_posterior(paramfile,post=post,$
;        isedfit_sfhgrid_dir=isedfit_sfhgrid_dir,iopath=isedpath,$
;        age=age,Z=Z,tau=tau,sfr0=sfr0,b100=b100,av=av,sfrage=sfrage)
;   endif

; ---------------------------------------------------------------------------
; posterior distributions: mass, sfrage, sfr100, av, Z, sSFR
    psfile = datapath+'santorini_posteriors.eps'
    im_plotconfig, 14, pos, psfile=psfile, height=5.0, ymargin=[1.0,1.0]

    plot, [0], [0], xsty=5, ysty=5, /noerase, /nodata, position=pos2[*,0], yrange=[0,1.05], $
      xrange=[6.2,9.8], ytitle='', charsize=1.4, xtitle=''
    im_plothist, mstar[*,ii]-alog10(cat[ii].mu), bin=0.15, /peak, /overplot, $
      /fill, fcolor=im_color('grey80')
    djs_oplot, (ised[ii].mass_50-alog10(cat[ii].mu))*[1,1], !y.crange, $
      line=5, thick=8
    plot, [0], [0], xsty=1, ysty=1, /noerase, /nodata, position=pos2[*,0], yrange=[0,1.05], $
      xrange=[6.2,9.8], ytitle='Probability', charsize=1.4, $
      xtitle=textoidl('log (M_{*}/M'+sunsymbol()+')'), $
      xtickinterval=1, ytickname=replicate(' ',10)

    xr = [-0.04,1] & bin = 0.05
;      xr = [-0.04,0.5] & bin = 0.02
    plot, [0], [0], xsty=5, ysty=5, /noerase, /nodata, yrange=[0,1.05], $
      xrange=xr, position=pos2[*,1], ytitle='', charsize=1.4, $
      xtitle=''
    im_plothist, sfrage[*,ii], bin=bin, /peak, /overplot, $
      /fill, fcolor=im_color('grey80')
    djs_oplot, median(sfrage[*,ii])*[1,1], !y.crange, line=5, thick=8
    plot, [0], [0], xsty=1, ysty=1, /noerase, /nodata, yrange=[0,1.05], $
      xrange=xr, position=pos2[*,1], ytitle='', charsize=1.4, $
      xtitle=textoidl('Age_{w} (Gyr)'), ytickname=replicate(' ',10), $
      xtickinterval=0.2

    plot, [0], [0], xsty=5, ysty=5, /noerase, /nodata, yrange=[0,1.1], $
      xrange=[-0.04,1.6], position=pos2[*,2], ytitle='', $
      xtitle='', charsize=1.4
    im_plothist, Z[*,ii]/0.02, bin=0.1, /peak, /overplot, $
      /fill, fcolor=im_color('grey80')
    djs_oplot, ised[ii].Z_50/0.02*[1,1], !y.crange, line=5, thick=8
    plot, [0], [0], xsty=1, ysty=1, /noerase, /nodata, yrange=[0,1.1], $
      xrange=[-0.04,1.6], position=pos2[*,2], ytitle='', $
      xtitle='Z/Z'+sunsymbol(), ytickname=replicate(' ',10), charsize=1.4, $
      xtickinterval=0.5

    xr = [-0.03,5.05] & yr = [0,1.1] & xtic = 1
;      xr = [-0.03,2.05] & yr = [0,1.1] & xtic = 0.5
    plot, [0], [0], xsty=5, ysty=5, /noerase, /nodata, yrange=yr, $
      xrange=xr, position=pos2[*,3], ytitle='', $
      xtitle=textoidl('A_{V} (mag)'), $
      ytickname=replicate(' ',10), charsize=1.4, xtickinterval=xtic
    im_plothist, av[*,ii], bin=0.3, /peak, /overplot, $
      /fill, fcolor=im_color('grey80')
    djs_oplot, ised[ii].av_50*[1,1], !y.crange, line=5, thick=8
    plot, [0], [0], xsty=1, ysty=1, /noerase, /nodata, yrange=yr, $
      xrange=xr, position=pos2[*,3], ytitle='Probability', $
      xtitle=textoidl('A_{V} (mag)'), $
      ytickname=replicate(' ',10), charsize=1.4, xtickinterval=xtic

    plot, [0], [0], xsty=5, ysty=5, /noerase, /nodata, yrange=[0,1.1], $
      xrange=[-0.5,2.5], position=pos2[*,4], ytitle='', $
      xtitle=textoidl('SFR (M'+sunsymbol()+' yr^{-1})'), $
      ytickname=replicate(' ',10), charsize=1.4, xtickinterval=2
    im_plothist, 10^(sfr0[*,ii]-alog10(cat[ii].mu)), bin=0.2, /peak, /overplot, $
      /fill, fcolor=im_color('grey80')
    djs_oplot, 10^(ised[ii].sfr_50-alog10(cat[ii].mu))*[1,1], !y.crange, line=5, thick=8
    plot, [0], [0], xsty=1, ysty=1, /noerase, /nodata, yrange=[0,1.1], $
      xrange=[-0.5,2.5], position=pos2[*,4], ytitle='', $
      xtitle=textoidl('SFR (M'+sunsymbol()+' yr^{-1})'), $
      ytickname=replicate(' ',10), charsize=1.4, xtickinterval=2

    plot, [0], [0], xsty=5, ysty=5, /noerase, /nodata, yrange=[0,1.1], $
      xrange=[0.0,10.0], position=pos2[*,5], ytitle='', $
      ytickname=replicate(' ',10), charsize=1.4
    im_plothist, 10^ssfr[*,ii], bin=0.6, /peak, /overplot, $
      /fill, fcolor=im_color('grey80')
    djs_oplot, median(10^ssfr[*,ii])*[1,1], !y.crange, line=5, thick=8
    plot, [0], [0], xsty=1, ysty=1, /noerase, /nodata, yrange=[0,1.1], $
      xrange=[0.0,10.0], position=pos2[*,5], ytitle='', $
      xtitle=textoidl('sSFR (Gyr^{-1})'), $
      ytickname=replicate(' ',10), charsize=1.4, xtickinterval=2

    im_plotconfig, /psclose, psfile=psfile, /pdf, /pskeep

stop

; ---------------------------------------------------------------------------
; SED plot for the paper
    paramfile = isedpath+'santorini_supergrid01_isedfit.par'
    model = isedfit_restore(paramfile,ised,iopath=isedpath,$
      isedfit_sfhgrid_dir=isedfit_sfhgrid_dir)

    psfile = datapath+'santorini_sed.eps'
    im_plotconfig, 0, pos, psfile=psfile, height=5.0, ymargin=[1.0,1.0]

    xtitle1 = textoidl('Observed-Frame Wavelength (\AA)')
    xtitle2 = textoidl('Rest-Frame Wavelength (\AA)')
    ytitle1 = textoidl('Apparent AB Magnitude')

    yrange = [31.0,24.0]
    xrange1 = [1900,60000]

    ticks1 = [4000,12000,30000]
    ticks2 = [250,500,1000,2000,4000]

    xrange2 = xrange1/(1.0+ised[0].zobj)
    plot, [0], [0], /nodata, xrange=xrange1, yrange=yrange, $
      xsty=9, ysty=1, xtitle=xtitle1, ytitle=ytitle1, xlog=1, $
      position=pos, xtickformat='(I0)', xticks=n_elements(ticks1)-1, xtickv=ticks1
    axis, /xaxis, xsty=1, xtitle=xtitle2, xrange=xrange2, xlog=1, $
      xticks=n_elements(ticks2)-1, xtickv=ticks2
;   im_legend, [cat[ii].galaxy,'z_{phot} = '+string(cat[0].z,format='(F4.2)')], $
;     /left, /top, box=0, margin=0, charsize=1.7
    oplot, model[0].wave, model[0].flux, line=0
    oplot, model[1].wave, model[1].flux, line=0, color=im_color('grey70')
    
    used = where((ised[0].maggies gt 0.0) and $ ; used in the fitting
      (ised[0].ivarmaggies gt 0.0),nused)
    upper = where((ised[0].maggies le 0.0) and $ ; upper limit
      (ised[0].ivarmaggies gt 0.0),nupper)

    if (nused ne 0L) then begin
       mab = maggies2mag(ised[0].maggies[used],$
         ivar=ised[0].ivarmaggies[used],magerr=mab_err)
       oploterror, weff[used], mab, mab_err, psym=symcat(16), $
         symsize=2.0, color=im_color('firebrick'), $
         errcolor=im_color('firebrick'), errthick=!p.thick
    endif

    if (nupper ne 0) then begin
       mab = maggies2mag(1.0/sqrt(ised[0].ivarmaggies[upper]))
       oploterror, weff[upper], mab, mab*0.0, psym=symcat(11,thick=6), $
         symsize=3.0, color=im_color('steel blue'), $
         errcolor=im_color('steel blue'), errthick=!p.thick
    endif

    im_plotconfig, /psclose, psfile=psfile, /pdf, /pskeep

stop    
    

return
end
