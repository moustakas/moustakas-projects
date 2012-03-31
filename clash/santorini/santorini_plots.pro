;;+
;; Calls ja_plot of ja_oplot but designed to do histograms
;; properly (not like psym=10. Especially useful with log axes
;; N.B. will loose top bin - set max to max(xval)+binsize to  avoid this
;; INPUTS:
;;  x - locations from histogram (i.e. the lower limit of the bin)
;;  y - histogram values
;;-
;pro ja_histplot,x,y,over=over,_extra=_extra
;
;  xr=reform(transpose(rebin(x,n_elements(x),2)),n_elements(x)*2)
;  yr=shift(reform(transpose(rebin(y,n_elements(y),2)),n_elements(y)*2),+1)
;  yr[0]=0
;  yr[n_elements(yr)-1]=0
;
;  if keyword_set(over) then ja_oplot,xr,yr,_extra=_extra $
;  else ja_plot,xr,yr,_extra=_extra
;end
;
;IDL> h=histogram(alog10(test),locations=l,binsize=0.5,max=4)
;IDL> ja_histplot,10^l,h,/xlog,yrange=[0,2]

pro render_postplot, xx, pos, xrange=xrange, yrange=yrange, $
  xtitle=xtitle, ytitle=ytitle, binsize=binsize, noerase=noerase, $
  xtickinterval=xtickinterval, nomedian=nomedian, monte=monte, $
  charsize=charsize, nonorm=nonorm, color_fill=color_fill, $
  color_outline=color_outline, color_monte=color_monte, $
  fill_monte=fill_monte, xlog=xlog, logbins=logbins, _extra=extra

    yrange = [0,1.1]
;   yrange = [0,1.05]

    if (n_elements(xrange) eq 0) then xrange = minmax(xx)*[0.9,1.1]
    if (n_elements(binsize) eq 0) then begin
       if keyword_set(logbins) then $
         binsize = (alog10(xrange[1])-alog10(xrange[0]))/ceil(0.3*sqrt(n_elements(xx))) else $
           binsize = (xrange[1]-xrange[0])/ceil(0.3*sqrt(n_elements(xx)))
    endif

    if n_elements(color_fill) eq 0 then color_fill = 'pale turquoise'
    if n_elements(color_outline) eq 0 then color_outline = 'pale turquoise'
    if n_elements(color_monte) eq 0 then color_monte = 'grey80'

    plot, [0], [0], xsty=5, ysty=5, /nodata, position=pos, $
      yrange=yrange, xrange=xrange, noerase=noerase, charsize=charsize, $
      xlog=xlog
    im_plothist, xx, bin=binsize, /peak, /overplot, /fill, $
      fcolor=im_color(color_fill), xhist, yhist, charsize=charsize, $
      color=im_color(color_outline,255), logbins=logbins, xlog=xlog
;   if keyword_set(nomedian) eq 0 then $
;     djs_oplot, median(xx)*[1,1], !y.crange, line=5, thick=8
    if n_elements(monte) ne 0 then begin
       im_plothist, monte, bin=binsize*2, mxhist, myhist, /noplot, $
         _extra=extra, logbins=logbins, xlog=xlog
       if keyword_set(nonorm) then normfactor = 1D else $
         normfactor = max(myhist)/(yrange[1]*0.9)
       if keyword_set(fill_monte) then begin
          im_plothist, monte, bin=binsize*2, /overplot, /fill, $
            normfactor=normfactor, fcolor=im_color(color_monte), $
            color=im_color(color_monte), logbins=logbins, xlog=xlog
       endif else begin
          im_plothist, monte, bin=binsize*2, /overplot, line=1, $
            normfactor=normfactor, color=im_color(color_monte), $
            logbins=logbins, xlog=xlog
;           normfactor=max(myhist)/max(yhist)/1.2
       endelse
    endif
    plot, [0], [0], xsty=1, ysty=1, /noerase, /nodata, position=pos, $
      yrange=yrange, xrange=xrange, ytitle=ytitle, $
      xtitle=textoidl(xtitle), xtickinterval=xtickinterval, $
      ytickname=replicate(' ',10), charsize=charsize, xlog=xlog, $
      _extra=extra
return
end

pro santorini_rebuild_posteriors

    common com_santorini, allpost, allised, postmodel, lowz_postmodel
    
    datapath = clash_path(/santorini)
    isedpath = datapath+'isedfit/'

    isedfit_sfhgrid_dir = datapath+'montegrids/'
    sfhgrid_paramfile = getenv('CLASH_DIR')+'/santorini/santorini_sfhgrid.par'
    params = yanny_readone(sfhgrid_paramfile)
    nmodel = params.nmonte
    ndraw = isedfit_ndraw()

    supergrid = [1,2,3,5,6] ; leave off the bursts
    super = get_santorini_supergrid(supergrid,nsuper=nsuper)
;   super = super[0:2] & nsuper = 3
    struct_print, super
    
    allpost = {supergrid: -1, mstar: fltarr(ndraw), sfrage: fltarr(ndraw), Z: fltarr(ndraw), $
      age: fltarr(ndraw), tau: fltarr(ndraw), sfr0: fltarr(ndraw), $
      av: fltarr(ndraw), chi2: fltarr(ndraw)};, bigsfr: fltarr(nmodel), bigsfrage: fltarr(nmodel), $
;     bigmass: fltarr(nmodel)}
    allpost = replicate(allpost,nsuper)
    allpost.supergrid = supergrid

    for ii = 0, nsuper-1 do begin
       if supergrid[ii] eq 5 then suffix = '_lowz' else suffix = ''

       paramfile = isedpath+'santorini'+suffix+'_supergrid0'+strtrim(supergrid[ii],2)+'_isedfit.par'
       splog, paramfile

       model = isedfit_restore(paramfile,allised1,iopath=isedpath,$
         isedfit_sfhgrid_dir=isedfit_sfhgrid_dir,/nomodel)
       if allised1.nburst eq 0 then begin
          if ii eq 0 then allised = allised1 else allised = [allised,allised1]
       endif

       delvarx, post
       mstar = isedfit_reconstruct_posterior(paramfile,post=post,$
         isedfit_sfhgrid_dir=isedfit_sfhgrid_dir,iopath=isedpath,$
         age=age,Z=Z,tau=tau,sfr0=sfr0,b100=b100,av=av,sfrage=sfrage,$
         chunkindx=chunkindx,modelindx=modelindx,indxage=ageindx,$
         bigsfr0=bigsfr,bigmass=bigmass,bigsfrage=bigsfrage)

       allpost[ii].mstar = mstar
       allpost[ii].sfrage = sfrage
       allpost[ii].Z = Z
       allpost[ii].age = age
       allpost[ii].tau = tau
       allpost[ii].sfr0 = sfr0
       allpost[ii].av = av
       allpost[ii].chi2 = post.chi2
       
;      allpost[ii].bigsfr = alog10(bigsfr)
;      allpost[ii].bigmass = alog10(bigmass)
;      allpost[ii].bigsfrage = bigsfrage

; reconstruct the posterior models, but just for the fiducial grid 
       if supergrid[ii] eq 6 then begin
          temp = replicate({zobj: allised1.zobj, chi2: 0.0, chunkindx: 0L, $
            modelindx: 0L, ageindx: 0L, scale: 0.0},ndraw)
          temp.chunkindx = chunkindx
          temp.modelindx = modelindx
          temp.ageindx = ageindx
          temp.scale = post.scale
;         temp.chi2 = post.chi2

          splog, paramfile
          postmodel = isedfit_restore(paramfile,in_isedfit=temp,$;/fnu,$
            isedfit_sfhgrid_dir=isedfit_sfhgrid_dir,iopath=isedpath)
;         plot, postmodel[0].wave, postmodel[0].flux, /xlog, xr=[1E4,4E4], yr=[30,20], ps=3
;         for bb = 1, ndraw-1 do djs_oplot, postmodel[bb].wave, postmodel[bb].flux, ps=3
       endif

; reconstruct the posterior models for the lowz grid
       if supergrid[ii] eq 5 then begin
          temp = replicate({zobj: allised1.zobj, chi2: 0.0, chunkindx: 0L, $
            modelindx: 0L, ageindx: 0L, scale: 0.0},ndraw)
          temp.chunkindx = chunkindx
          temp.modelindx = modelindx
          temp.ageindx = ageindx
          temp.scale = post.scale
;         temp.chi2 = post.chi2

          splog, paramfile
          lowz_postmodel = isedfit_restore(paramfile,in_isedfit=temp,$;/fnu,$
            isedfit_sfhgrid_dir=isedfit_sfhgrid_dir,iopath=isedpath)
;         plot, postmodel[0].wave, postmodel[0].flux, /xlog, xr=[1E4,4E4], yr=[30,20], ps=3
;         for bb = 1, ndraw-1 do djs_oplot, postmodel[bb].wave, postmodel[bb].flux, ps=3
       endif
    endfor

return
end    

pro santorini_plots, supergrid, models=models, isedfit=isedfit, $
  qaplot=qaplot, clobber=clobber
; jm11nov08ucsd - build plots for the paper

    common com_santorini

    datapath = clash_path(/santorini)
    isedpath = datapath+'isedfit/'

    isedfit_sfhgrid_dir = datapath+'montegrids/'
    sfhgrid_paramfile = getenv('CLASH_DIR')+'/santorini/santorini_sfhgrid.par'

    filters = clash_filterlist(nice=nice_filters,/useirac,pivotwave=pivotwave,width=width)
    ndraw = isedfit_ndraw()

; gather the photometry and restore the results
    cat = read_santorini()
    nobj = n_elements(cat)
    logmu = alog10(cat.mu)
    
    if n_elements(allpost) eq 0 then santorini_rebuild_posteriors

    jj = mrdfits(isedpath+'santorini_lowz_fsps_chab_calzetti_sfhgrid05.fits.gz',1)
    niceprint, allpost.supergrid, 10^(allised.mass_50-logmu), 10^(allised.sfr_50-logmu), $
      allised.sfrage_50*1E3, 10^(allised.mass-logmu), 10^(allised.sfr-logmu), $
      allised.sfrage*1E3, allised.tau*1E3
;   niceprint, [allised.chi2,jj.chi2]

    splog, 'Supergrid 1 - free dust, metals'
    tsfr = weighted_quantile(allpost[0].sfrage*1E3,quant=0.95)
    splog, tsfr, 1D3*allised[0].sfrage, getredshift(getage(9.6)-tsfr/1E3), $
      getredshift(getage(9.6)-allised[0].sfrage) 
    
    splog, 'Supergrid 6 - fiducial'
    tsfr = weighted_quantile(allpost[4].sfrage*1E3,quant=0.95)
    splog, tsfr, 1D3*allised[4].sfrage, getredshift(getage(9.6)-tsfr/1E3), $
      getredshift(getage(9.6)-allised[4].sfrage)
;   splog, tsfr, lf_t2z(0.7*(lf_z2t(9.6)/0.7-tsfr/1D3)

    splog, 'Supergrid 3 - simple tau'
    tsfr = weighted_quantile(allpost[2].sfrage*1E3,quant=0.95)
    splog, tsfr, 1D3*allised[2].sfrage, getredshift(getage(9.6)-tsfr/1E3), $
      getredshift(getage(9.6)-allised[2].sfrage)

; ---------------------------------------------------------------------------
; SED plot for the paper
    paramfile = isedpath+'santorini_supergrid06_isedfit.par'
    model = isedfit_restore(paramfile,ised,iopath=isedpath,$
      isedfit_sfhgrid_dir=isedfit_sfhgrid_dir);,/fnu)

    psfile = datapath+'santorini_sed.eps'
    im_plotconfig, 0, pos, psfile=psfile, height=5.0, ymargin=[1.0,1.0]

    xtitle1 = textoidl('Observed-frame wavelength (\mu'+'m)')
    xtitle2 = textoidl('Rest-frame wavelength (\mu'+'m)')
;   ytitle1 = textoidl('Flux Density (10^{-30} '+fnu_units()+')')
    ytitle1 = textoidl('AB magnitude')

;   scale = 1D30
;   yrange = [-0.5,6]
    yrange = [29.5,23.0]
    xrange1 = [0.19,7.0]

;   ticks1 = [0.3,1.2,4]
    ticks1 = [0.2,0.3,0.4,0.6,0.8,1,2,3,4,5,6]
    ticks2 = [0.02,0.04,0.06,0.1,0.2,0.3,0.4,0.6]

    xrange2 = xrange1/(1.0+ised.zobj)
    plot, [0], [0], /nodata, xrange=xrange1, yrange=yrange, $
      xsty=5, ysty=5, xlog=1, position=pos
    polyfill, [postmodel.wave/1D4,reverse(postmodel.wave/1D4)], $
      [min(postmodel.flux,dim=2),reverse(max(postmodel.flux,dim=2))], $
      /data, color=im_color('pale turquoise'), noclip=0, /fill
    plot, [0], [0], /nodata, /noerase, xrange=xrange1, yrange=yrange, $
      xsty=9, ysty=1, xtitle=xtitle1, ytitle=ytitle1, xlog=1, $
      position=pos, xticks=n_elements(ticks1)-1, xtickv=ticks1, $
      xtickformat='(G0)', yminor=2
    
;    plot, [0], [0], /nodata, xrange=xrange1, yrange=yrange, $
;      xsty=9, ysty=1, xtitle=xtitle1, ytitle=ytitle1, xlog=1, $
;      position=pos, xtickformat='(I0)', xticks=n_elements(ticks1)-1, xtickv=ticks1
;    nww = 300
;    ww = shuffle_indx(ndraw,num=nww)
;    for ii = 0, nww-1 do begin
;       djs_oplot, postmodel[ww[ii]].wave/1D4, postmodel[ww[ii]].flux, $
;         line=0, thick=1, color=im_color('black')
;;      cc = get_kbrd(1)
;    endfor

    axis, /xaxis, xsty=1, xtitle='', xrange=xrange2, xlog=1, $
      xticks=n_elements(ticks2)-1, xtickv=ticks2, xtickformat='(G0)'
    xyouts, (pos[2]-pos[0])/2+pos[0], pos[3]+0.07, textoidl(xtitle2), $
      align=0.5, /norm
    oplot, model.wave/1D4, model.flux, line=0, thick=5;, color=im_color('brown4')

;    snr = ised.maggies*sqrt(ised.ivarmaggies) ; fractional error
;    good = where(snr ge 2.0,comp=upper)
;
;    bestflux = scale*ised.bestmaggies*10^(-0.4*48.6)
;    flux = scale*ised.maggies[good]*10^(-0.4*48.6)
;    ferr = scale*10^(-0.4*48.6)/sqrt(ised.ivarmaggies[good])
;    fupper = scale*2.0/sqrt(ised.ivarmaggies[upper])*10^(-0.4*48.6)
;    
;    djs_oplot, pivotwave, bestflux, psym=symcat(6,thick=6), $
;      symsize=2.5, color=im_color('steel blue')
;    oploterror, pivotwave[good], flux, ferr, psym=symcat(16), $
;      symsize=2.0, color=im_color('firebrick'), $
;      errcolor=im_color('firebrick'), errthick=8
;    djs_oplot, pivotwave[upper], fupper, psym=symcat(11,thick=6), $
;      symsize=3.0, color=im_color('grey20')
       
    mab = maggies2mag(ised.maggies,ivar=ised.ivarmaggies,$
      lomagerr=loerr,himagerr=hierr,magnsigma=mabupper,nsigma=1.0)
    splog, 'Hack!'
    mabupper[1] = mab[1] & mab[1] = -99 
    good = where(mab gt -99.0,ngood,comp=upper,ncomp=nupper)
    
    djs_oplot, pivotwave/1D4, -2.5*alog10(ised.bestmaggies), $
      psym=symcat(6,thick=6), symsize=3, color=im_color('firebrick4')

    col = 'navy'
    oploterror, pivotwave[good]/1D4, mab[good], width[good]/1D4/2, loerr[good], psym=symcat(16), $
      symsize=2.2, color=im_color(col), $
      errcolor=im_color(col), errthick=8, /lobar
    oploterror, pivotwave[good]/1D4, mab[good], width[good]/1D4/2, hierr[good], psym=3, $
      color=im_color(col), errcolor=im_color(col), errthick=8, /hibar
;   oploterror, pivotwave[upper]/1D4, mabupper[upper], width[upper]/1D4/2, mabupper[upper]*0, $
;     psym=symcat(11,thick=6), symsize=2.4, color=im_color('dark green'), $
;     errcolor=im_color('dark green')
    djs_oplot, pivotwave[upper]/1D4, mabupper[upper], psym=symcat(11,thick=6), $
      symsize=2.4, color=im_color('dark green')

; show an inset with the posterior distribution on <t>_SFR
    showinset = 0
    
    if showinset then begin
       csize = 1.4
       pos1 = [0.23,0.6,0.52,0.82]
       render_postplot, allpost[1].sfrage*1E3, pos1, xrange=[-10,400], $
         /noerase, ytitle='', xtitle='<t>_{SFR} (Myr)', $
         xtickinterval=100, charsize=csize
       xyouts, pos1[0]-0.02, (pos1[3]-pos1[1])/2+pos1[1], 'Probability', $
         align=0.5, orientation=90, charsize=csize, /norm
    endif
    
    im_plotconfig, /psclose, psfile=psfile, /pdf, /pskeep

; ---------------------------------------------------------------------------
; SED plot for the paper - low-redshift solution
    paramfile = isedpath+'santorini_lowz_supergrid05_isedfit.par'
    model = isedfit_restore(paramfile,ised,iopath=isedpath,$
      isedfit_sfhgrid_dir=isedfit_sfhgrid_dir);,/fnu)

    psfile = datapath+'santorini_lowz_sed.eps'
    im_plotconfig, 0, pos, psfile=psfile, height=5.0, ymargin=[1.0,1.0]

    xtitle1 = textoidl('Observed-frame wavelength (\mu'+'m)')
    xtitle2 = textoidl('Rest-frame wavelength (\mu'+'m)')
;   ytitle1 = textoidl('Flux Density (10^{-30} '+fnu_units()+')')
    ytitle1 = textoidl('AB magnitude')

    scale = 1D30
    yrange = [29.5,23.0]
;   yrange = [-0.5,9]
    xrange1 = [0.19,7.0]

    ticks1 = [0.2,0.3,0.4,0.6,0.8,1,2,3,4,5,6]
;   ticks2 = [0.02,0.04,0.06,0.1,0.2,0.3,0.4,0.6]
    ticks2 = [0.05,0.07,0.1,0.2,0.3,0.5,0.8,1]

    xrange2 = xrange1/(1.0+ised.zobj)
    plot, [0], [0], /nodata, xrange=xrange1, yrange=yrange, $
      xsty=5, ysty=5, xlog=1, position=pos
    polyfill, [lowz_postmodel.wave/1D4,reverse(lowz_postmodel.wave/1D4)], $
      [min(lowz_postmodel.flux,dim=2),reverse(max(lowz_postmodel.flux,dim=2))], $
      /data, color=im_color('grey80'), noclip=0, /fill
    plot, [0], [0], /nodata, /noerase, xrange=xrange1, yrange=yrange, $
      xsty=9, ysty=1, xtitle=xtitle1, ytitle=ytitle1, xlog=1, $
      position=pos, xticks=n_elements(ticks1)-1, xtickv=ticks1, $
      xtickformat='(G0)', yminor=2

;    plot, [0], [0], /nodata, xrange=xrange1, yrange=yrange, $
;      xsty=9, ysty=1, xtitle=xtitle1, ytitle=ytitle1, xlog=1, $
;      position=pos, xtickformat='(G0)', xticks=n_elements(ticks1)-1, xtickv=ticks1

    axis, /xaxis, xsty=1, xtitle='', xrange=xrange2, xlog=1, $
      xticks=n_elements(ticks2)-1, xtickv=ticks2, xtickformat='(G0)'
    xyouts, (pos[2]-pos[0])/2+pos[0], pos[3]+0.07, textoidl(xtitle2), $
      align=0.5, /norm
    oplot, model.wave/1D4, model.flux, line=0, thick=5
;   oplot, model.wave, scale*model.flux, line=0

;    snr = ised.maggies*sqrt(ised.ivarmaggies) ; fractional error
;    good = where(snr ge 1.0,comp=upper)
;    bestflux = scale*ised.bestmaggies*10^(-0.4*48.6)
;;   bestflux = scale*ised.bestmaggies*10^(-0.4*48.6)
;    flux = scale*ised.maggies[good]*10^(-0.4*48.6)
;    ferr = scale*10^(-0.4*48.6)/sqrt(ised.ivarmaggies[good])
;    fupper = scale*2.0/sqrt(ised.ivarmaggies[upper])*10^(-0.4*48.6)
;    
;    djs_oplot, pivotwave, bestflux, psym=symcat(6,thick=6), $
;      symsize=2.5, color=im_color('steel blue')
;    oploterror, pivotwave[good], flux, ferr, psym=symcat(16), $
;      symsize=2.0, color=im_color('firebrick'), $
;      errcolor=im_color('firebrick'), errthick=8
;    djs_oplot, pivotwave[upper], fupper, psym=symcat(11,thick=6), $
;      symsize=3.0, color=im_color('grey20')

    mab = maggies2mag(ised.maggies,ivar=ised.ivarmaggies,$
      lomagerr=loerr,himagerr=hierr,magnsigma=mabupper,nsigma=1.0)
    splog, 'Hack!'
    mabupper[1] = mab[1] & mab[1] = -99 
    good = where(mab gt -99.0,ngood,comp=upper,ncomp=nupper)

    djs_oplot, pivotwave/1D4, -2.5*alog10(ised.bestmaggies), $
      psym=symcat(6,thick=6), symsize=2.5, color=im_color('firebrick')

    col = 'navy'
    oploterror, pivotwave[good]/1D4, mab[good], width[good]/1D4/2, loerr[good], psym=symcat(16), $
      symsize=2.2, color=im_color(col), $
      errcolor=im_color(col), errthick=8, /lobar
    oploterror, pivotwave[good]/1D4, mab[good], width[good]/1D4/2, hierr[good], psym=3, $
      color=im_color(col), errcolor=im_color(col), errthick=8, /hibar

    djs_oplot, pivotwave[upper]/1D4, mabupper[upper], psym=symcat(11,thick=6), $
      symsize=2.4, color=im_color('dark green')
    
; show an inset with the posterior probabilities of the low- and
; high-z solutions
    csize = 1.4
    pos1 = [0.23,0.6,0.52,0.82]
    bin = 0.5
    
    render_postplot, allpost[1].chi2, pos1, xrange=[0,25], $
      /noerase, ytitle='', xtitle='\chi^{2}', $
      charsize=csize, binsize=bin, monte=allpost[3].chi2, $
      color_monte='grey80', color_fill='pale turquoise', $
      /fill_monte
;   im_plothist, , binsize=bin, /over
    xyouts, pos1[0]-0.02, (pos1[3]-pos1[1])/2+pos1[1], 'Frequency', $
      align=0.5, orientation=90, charsize=csize, /norm
    
    im_plotconfig, /psclose, psfile=psfile, /pdf, /pskeep

; ---------------------------------------------------------------------------
; posterior distributions plots - supergrid06
    monte = mrdfits(isedfit_sfhgrid_dir+'sfhgrid06/fsps/chab_montegrid.fits.gz',1)

    psfile = datapath+'santorini_posteriors.eps'
    im_plotconfig, 3, pos, psfile=psfile, xspace=0.2*[1,1], yspace=0.0, $
      xmargin=[0.8,0.3], height=3.0, width=3*[1,1,1], charsize=1.6, $
      color_outline='black'

    render_postplot, 10^(allpost[1].mstar-logmu), pos[*,0], xrange=[1D7,3D9], $
      ytitle='', /logbins, binsize=0.1, /xlog, $
      xtitle='Stellar Mass [(\mu/15)^{-1} M'+sunsymbol()+']', $
      color_outline='black'
    xyouts, pos[0,0]-0.02, (pos[3,0]-pos[1,0])/2.0+pos[1,0], 'Probability', $
      align=0.5, orientation=90, /norm
    
    render_postplot, 10^(allpost[1].sfr0-logmu), pos[*,1], /noerase, xrange=[-0.1,2.3], $
      xtitle='SFR [(\mu/15)^{-1} M'+sunsymbol()+' yr^{-1}]', bin=0.1, $
      color_outline='black' ;, monte=allpost[1].bigsfr+10
    
    render_postplot, allpost[1].sfrage*1E3, pos[*,2], xrange=[-10,400], /noerase, $
      xtitle='<t>_{SFR} (Myr)', xtickinterval=100, color_outline='black', bin=20

    im_plotconfig, /psclose, psfile=psfile, /pdf, /pskeep


stop
stop
    
; ---------------------------------------------------------------------------
; SFRAGE vs AV for supergrid01 for Wei
    psfile = datapath+'santorini_sfrage_vs_av.eps'
    im_plotconfig, 0, pos, psfile=psfile, height=5.0, ymargin=[1.0,1.0]

    hogg_scatterplot, allpost[0].sfrage*1E3, allpost[0].av, $
      xrange=[-20,300], yrange=[-0.1,2], xsty=1, ysty=1, $
      xtitle=textoidl('<t>_{SFR} (Myr)'), ytitle=textoidl('A_{V} (mag)'), $
      /outliers, levels=[0.25,0.5,0.75,0.9], outcolor=djs_icolor('black'), $
      ynpix=12, xnpix=12, /nogrey, /internal
    
    im_plotconfig, /psclose, psfile=psfile, /pdf, /pskeep    

; ---------------------------------------------------------------------------
; compare the posterior distributions on mass, SFR, and age
    psfile = datapath+'santorini_sfh_dependence.eps'
    im_plotconfig, 5, pos, psfile=psfile, xspace=0.3, yspace=1.0

    c1 = 'black'
    c2 = 'firebrick'
    c3 = 'dodger blue'

    l1 = 0
    l2 = 5
    l3 = 3
    
; mass
    djs_plot, [0], [0], xsty=1, ysty=1, /nodata, yrange=[0,1.1], $
      xrange=[7.0,10.0], position=pos[*,0], ytitle='Probability', $
      xtitle='log (M_{*}/M'+sunsymbol()+')', xtickinterval=1, $
      ytickname=replicate(' ',10)
    bin = (!x.crange[1]-!x.crange[0])/ceil(0.3*sqrt(ndraw))
    im_plothist, allpost[1].mstar-logmu, bin=bin, /peak, /overplot, color=im_color(c1), thick=8, line=l1
    im_plothist, allpost[2].mstar-logmu, bin=bin, /peak, /overplot, color=im_color(c2), thick=6, line=l2
    im_plothist, allpost[3].mstar-logmu, bin=bin, /peak, /overplot, color=im_color(c3), thick=6, line=l3

; SFR    
    djs_plot, [0], [0], xsty=1, ysty=1, /nodata, /noerase, yrange=[0,1.1], $
      xrange=[-1.2,1], position=pos[*,1], ytitle='', $
      xtitle='log (SFR /M'+sunsymbol()+' yr^{-1})', ytickname=replicate(' ',10)
    bin = (!x.crange[1]-!x.crange[0])/ceil(0.3*sqrt(ndraw))
    im_plothist, allpost[1].sfr0-logmu, bin=bin, /peak, /overplot, color=im_color(c1), thick=8, line=l1
    im_plothist, allpost[2].sfr0-logmu, bin=bin, /peak, /overplot, color=im_color(c2), thick=6, line=l2
    im_plothist, allpost[3].sfr0-logmu, bin=bin, /peak, /overplot, color=im_color(c3), thick=6, line=l3

; SFR-weighted AGE
    djs_plot, [0], [0], xsty=1, ysty=1, /nodata, /noerase, yrange=[0,1.1], $
      xrange=[-20,500], position=pos[*,2], ytitle='Probability', xtitle=textoidl('t_{w} (Myr)'), $
      xtickinterval=100, ytickname=replicate(' ',10)
    bin = (!x.crange[1]-!x.crange[0])/ceil(0.3*sqrt(ndraw))
    im_plothist, allpost[1].sfrage*1E3, bin=bin, /peak, /overplot, color=im_color(c1), thick=8, line=l1
    im_plothist, allpost[2].sfrage*1E3, bin=bin, /peak, /overplot, color=im_color(c2), thick=6, line=l2
    im_plothist, allpost[3].sfrage*1E3, bin=bin, /peak, /overplot, color=im_color(c3), thick=6, line=l3
    
    djs_plot, [0], [0], xsty=5, ysty=5, /nodata, /noerase, $
      position=pos[*,3], xrange=[0,1], yrange=[0,1]
    im_legend, ['Delayed \tau','Simple \tau','Bursty'], box=1, position=[0.15,0.7], $
      color=[c1,c2,c3], line=[l1,l2,l3], pspacing=1.9, thick=8

    im_plotconfig, /psclose, psfile=psfile, /pdf, /pskeep

; ---------------------------------------------------------------------------
; posterior distributions plots - supergrid01    
    monte = mrdfits(isedfit_sfhgrid_dir+'sfhgrid01/fsps/calzetti/chab_montegrid.fits.gz',1)

    psfile = datapath+'santorini_posteriors_dustymetals.eps'
    im_plotconfig, 17, pos, psfile=psfile, xspace=[0.3,0.3], yspace=1.0, $
      xmargin=[1.0,0.5], height=2.1*[1,1], width=2.7*[1,1,1], charsize=1.4

    render_postplot, allpost[0].mstar-logmu, pos[*,0], xrange=[7.0,10.0], $
      ytitle='Probability', xtitle='log (M_{*}/M'+sunsymbol()+')', $
      xtickinterval=1
;     monte=allpost[0].bigmass+median(allpost[0].mstar-logmu)+1.5
;   djs_oplot, !x.crange, !y.crange[1]*0.85*[1,1], line=1

    render_postplot, allpost[0].Z/0.019, pos[*,1], xrange=[-0.1,1.8], $
      xtitle='Z/Z'+sunsymbol(), /noerase, /nomedian ;, monte=monte.z/0.019
    
    render_postplot, allpost[0].av, pos[*,2], xrange=minmax(monte.av)+[-0.05,0.05], $
      xtitle='A_{V} (mag)', /noerase, xtickinterval=1;, monte=monte.av
    
    render_postplot, allpost[0].tau, pos[*,3], /noerase, xrange=[-0.2,1.2], $
      ytitle='Probability', xtitle='\tau (Gyr)', xtickinterval=0.5;monte=monte.tau, $
    
    render_postplot, allpost[0].sfr0-logmu, pos[*,4], /noerase, xrange=[-0.3,2.2], $ ; xrange=minmax(allpost[0].sfr0-logmu)*[0.99,1.01], $
      xtitle='log (SFR /M'+sunsymbol()+' yr^{-1})';, monte=allpost[0].bigsfr+10
    
    render_postplot, allpost[0].sfrage*1E3, pos[*,5], xrange=[-5,300], /noerase, $
      xtitle='t_{w} (Myr)', xtickinterval=100;, monte=allpost[0].bigsfrage*1E3

    im_plotconfig, /psclose, psfile=psfile, /pdf, /pskeep

;; ---------------------------------------------------------------------------
;; plot the NxN distributions of parameters as an upper triangle
;    this = 0 ; z=9.56 solution
;    
;    psfile = datapath+'santorini_manyd.eps'
;    splog, 'Writing '+psfile
;    manyd = transpose([[allpost[this].mstar],[allpost[this].sfrage],[allpost[this].Z/0.019],[allpost[this].av],[allpost[this].sfr0]])
;    label = textoidl(['log (M/M'+sunsymbol()+')','t_{w} (Gyr)','Z/Z'+sunsymbol(),'A_{V} (mag)','log (\psi/M'+sunsymbol()+' yr^{-1})'])
;    
;    im_manyd_scatterplot, fltarr(ndraw)+1, manyd, psfile, label=label, $
;      axis_char_scale=1.4, /internal, outliers=1, $
;      /nogrey, levels=errorf((dindgen(2)+1)/sqrt(2)), /upper, /fill
;    spawn, 'ps2pdf '+psfile+' '+repstr(psfile,'.eps','.pdf'), /sh

stop    
    
; ---------------------------------------------------------------------------
; posterior distributions plots - supergrid06
    monte = mrdfits(isedfit_sfhgrid_dir+'sfhgrid06/fsps/chab_montegrid.fits.gz',1)

    psfile = datapath+'santorini_posteriors.eps'
    im_plotconfig, 5, pos, psfile=psfile, xspace=0.3, yspace=1.0, $
      xmargin=[1.0,0.5], height=3.5*[1,1], width=4*[1,1], charsize=1.8, $
      color_outline='black'

    render_postplot, 10^(allpost[1].mstar-logmu), pos[*,0], xrange=[1D7,3D9], $
      ytitle='Probability', /logbins, binsize=0.1, /xlog, $
      xtitle='Stellar Mass [(\mu/15)^{-1} M'+sunsymbol()+']', $
      color_outline='black'

;    render_postplot, 10^(allpost[1].mstar-logmu), pos[*,0], xrange=10^[6.9,9.5], $
;      ytitle='Probability', /logbins, /xlog, binsize=0.1, $
;      xtitle='log_{10}(M_{*}/M'+sunsymbol()+')', $
;;     xtitle='log_{10}(M_{*}/M'+sunsymbol()+')-log_{10}(\mu/15)', $
;      color_outline='black'
;;   djs_oplot, (allised[1].mass-logmu)*[1,1], !y.crange, line=5, thick=8

    render_postplot, 10^(allpost[1].sfr0-logmu), pos[*,1], /noerase, xrange=[-0.1,2.5], $
;     xtitle='SFR (M'+sunsymbol()+' yr^{-1})', $
      xtitle='SFR [(\mu/15)^{-1} M'+sunsymbol()+' yr^{-1}]', $
      color_outline='black' ;, monte=allpost[1].bigsfr+10
;   djs_oplot, (allised[1].sfr-logmu)*[1,1], !y.crange, line=5, thick=8
    
    render_postplot, allpost[1].tau*1E3, pos[*,2], /noerase, xrange=[-100,1100], $
      ytitle='Probability', xtitle='\tau (Myr)', xtickinterval=500, $
      color_outline='black' ;monte=monte.tau, $
;   djs_oplot, allised[1].tau*1E3*[1,1], !y.crange, line=5, thick=8
    
    render_postplot, allpost[1].sfrage*1E3, pos[*,3], xrange=[-10,400], /noerase, $
      xtitle='<t>_{SFR} (Myr)', xtickinterval=100, color_outline='black'
;   im_plothist, allpost[1].sfrage*1E3, /cumu, /overplot, color='red'
;   djs_oplot, allised[1].sfrage*1E3*[1,1], !y.crange, line=5, thick=8

    im_plotconfig, /psclose, psfile=psfile, /pdf, /pskeep

;   outfile = datapath+'santorini_sed1.pdf'
;   allfiles = datapath+['santorini_supergrid02_sed.eps','santorini_supergrid02_posteriors.eps']
;   spawn, 'gs -q -dNOPAUSE -sDEVICE=pdfwrite '+$
;     '-dEPSFitPage -sOutputFile='+outfile+' -dBATCH '+strjoin(allfiles,' '), /sh

return
end
