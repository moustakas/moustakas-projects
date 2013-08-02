pro render_postplot, xx, pos, xrange=xrange, yrange=yrange, $
  xtitle=xtitle, ytitle=ytitle, binsize=binsize, noerase=noerase, $
  xtickinterval=xtickinterval, nomedian=nomedian, monte=monte, $
  charsize=charsize, nonorm=nonorm, color_fill=color_fill, $
  color_outline=color_outline, color_monte=color_monte, $
  fill_monte=fill_monte, xlog=xlog, logbins=logbins, _extra=extra, $
  overplot=overplot, nofill=nofill, linestyle=linestyle, thick=thick

    if n_elements(yrange) eq 0 then yrange = [0,1.1]
;   yrange = [0,1.05]

    if (n_elements(xrange) eq 0) then xrange = minmax(xx)*[0.9,1.1]
    if (n_elements(binsize) eq 0) then begin
       if keyword_set(logbins) then $
         binsize = (alog10(xrange[1])-alog10(xrange[0]))/ceil(0.3*sqrt(n_elements(xx))) else $
           binsize = (xrange[1]-xrange[0])/ceil(0.3*sqrt(n_elements(xx)))
    endif

    if n_elements(color_fill) eq 0 then color_fill = 'powder blue'
    if n_elements(color_outline) eq 0 then color_outline = 'powder blue'
    if n_elements(color_monte) eq 0 then color_monte = 'grey80'

    if keyword_set(overplot) eq 0 then begin
       plot, [0], [0], xsty=5, ysty=5, /nodata, position=pos, $
         yrange=yrange, xrange=xrange, noerase=noerase, charsize=charsize, $
         xlog=xlog
    endif
    im_plothist, xx, bin=binsize, peak=1, /overplot, fill=(keyword_set(nofill) eq 0), $
      fcolor=im_color(color_fill), xhist, yhist, charsize=charsize, $
      color=im_color(color_outline,255), logbins=logbins, xlog=xlog, $
      linestyle=linestyle, thick=thick
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
    if keyword_set(overplot) eq 0 then begin
       plot, [0], [0], xsty=1, ysty=1, /noerase, /nodata, position=pos, $
         yrange=yrange, xrange=xrange, ytitle=ytitle, $
         xtitle=textoidl(xtitle), xtickinterval=xtickinterval, $
         ytickname=replicate(' ',10), charsize=charsize, xlog=xlog, $
         _extra=extra
    endif
return
end

pro rebuild_posteriors

    common com_rebuild, allmodel_super1, allmodel_super2, allised_super1, $
      allised_super2, allpost_super1, allpost_super2, postmodel_super1, $
      postmodel_super2

    light = 2.99792458D18       ; speed of light [A/s]
    betawave = [15369.000,35510.001] ; pivot wavelengths for F160W,[ch1]

    prefix = 'z11'
    rootpath = getenv('IM_PROJECTS_DIR')+'/clash/z11/'
    isedfit_dir = rootpath
    montegrids_dir = rootpath+'montegrids/'
    isedfit_paramfile = isedfit_dir+prefix+'_paramfile.par'
    sfhgrid_paramfile = isedfit_dir+prefix+'_sfhgrid.par'
    supergrid_paramfile = isedfit_dir+prefix+'_supergrid.par'

    cat = read_z11()
    ngal = n_elements(cat)
    
    params = yanny_readone(sfhgrid_paramfile)
    nmodel = params.nmonte
    ndraw = isedfit_ndraw()

    super = read_supergrid_paramfile(supergrid_paramfile)
    struct_print, super

; restore the models and the best-fitting results
    allmodel_super1 = isedfit_restore(isedfit_paramfile,allised_super1,/fnu,$
      supergrid_paramfile=supergrid_paramfile,thissupergrid=1,$
      isedfit_dir=isedfit_dir,montegrids_dir=montegrids_dir)
    allmodel_super1.flux = allmodel_super1.flux*1D29 ; convert to uJy

    allmodel_super2 = isedfit_restore(isedfit_paramfile,allised_super2,/fnu,$
      supergrid_paramfile=supergrid_paramfile,thissupergrid=2,$
      isedfit_dir=isedfit_dir,montegrids_dir=montegrids_dir)
    allmodel_super2.flux = allmodel_super2.flux*1D29 ; convert to uJy

; build the posteriors    
    allpost = {mstar: fltarr(ndraw), sfrage: fltarr(ndraw), $
      sfrm: fltarr(ndraw), Z: fltarr(ndraw), age: fltarr(ndraw), tau: fltarr(ndraw), $
      sfr0: fltarr(ndraw), chi2: fltarr(ndraw), scale: fltarr(ndraw), $
      beta: fltarr(ndraw)}

; reconstruct the posterior models: supergrid01    
    allpost_super1 = replicate(allpost,ngal)
    mstar = isedfit_reconstruct_posterior(isedfit_paramfile,$
      post=post_super1,supergrid_paramfile=supergrid_paramfile,$
      thissupergrid=1,isedfit_dir=isedfit_dir,$
      chunkindx=chunkindx,modelindx=modelindx,indxage=ageindx,$
      montegrids_dir=montegrids_dir,age=age,Z=Z,tau=tau,$
      sfr0=sfr0,av=av,sfrage=sfrage)
    allpost_super1.mstar = mstar
    allpost_super1.sfrage = sfrage
    allpost_super1.age = age
    allpost_super1.sfrm = sfr0-mstar+9
    allpost_super1.sfr0 = sfr0
    allpost_super1.Z = Z
    allpost_super1.tau = tau
    allpost_super1.chi2 = post_super1.chi2
    allpost_super1.scale = post_super1.scale

    for ii = 0, ngal-1 do begin
       temp = replicate({zobj: allised_super1[0].zobj, chi2: 0.0, chunkindx: 0L, $
         modelindx: 0L, ageindx: 0L, scale: 0.0},ndraw)
       temp.chunkindx = chunkindx[*,ii]
       temp.modelindx = modelindx[*,ii]
       temp.ageindx = ageindx[*,ii]
       temp.scale = allpost_super1[ii].scale
          
       postmodel1 = isedfit_restore(isedfit_paramfile,in_isedfit=temp,/fnu,$
         supergrid_paramfile=supergrid_paramfile,thissupergrid=1,$
         isedfit_dir=isedfit_dir,montegrids_dir=montegrids_dir)
       postmodel1.flux = postmodel1.flux*1D29 ; convert to uJy
       if ii eq 0 then postmodel_super1 = postmodel1 else postmodel_super1 = [[postmodel_super1],[postmodel1]]

;      plot, postmodel[0].wave/1D4, postmodel[0].flux, /xlog, xr=[1,7], yr=[-0.05,0.7], ps=3
;      for bb = 1, ndraw-1 do djs_oplot, postmodel[bb].wave/1D4, postmodel[bb].flux, ps=3
          
; get the posterior for beta; for broadband filters beta can be
; derived from the effective wavelengths of the filters and the color
; (see, e.g., Dunlop+12, Section 3):
;
;   beta = F160W-[ch1]/(-2.5*alog10(F160W_wave/[ch2]_wave)) - 2

       coeff = -2.5*alog10(betawave[0]/betawave[1])
       for jj = 0L, ndraw-1 do begin
          flam = postmodel1[jj].flux*1D-29*light/postmodel1[jj].wave^2
          mm = reform(k_project_filters(k_lambda_to_edges(postmodel1[jj].wave),$
            flam,filterlist=['clash_wfc3_f160w.par','spitzer_irac_ch1.par'],$
            /silent))
;         mmlam = reform(mm)*10^(-0.4*48.6)*light/betawave^2
          color = -2.5*alog10(mm[0]/mm[1])
          beta1 = color/coeff-2.0

;         djs_plot, postmodel1[0].wave/1D4, flam, xr=[1,6], xsty=1
;         djs_oplot, betawave/1D4, mmlam, psym=8, color='green'
          allpost_super1[ii].beta[jj] = beta1
       endfor
    endfor 
    
; reconstruct the posterior models: supergrid02
    allpost_super2 = replicate(allpost,ngal)
    mstar = isedfit_reconstruct_posterior(isedfit_paramfile,$
      post=post_super2,supergrid_paramfile=supergrid_paramfile,$
      thissupergrid=2,isedfit_dir=isedfit_dir,$
      chunkindx=chunkindx,modelindx=modelindx,indxage=ageindx,$
      montegrids_dir=montegrids_dir,age=age,Z=Z,tau=tau,$
      sfr0=sfr0,av=av,sfrage=sfrage)
    allpost_super2.mstar = mstar
    allpost_super2.sfrage = sfrage
    allpost_super2.age = age
    allpost_super2.sfrm = sfr0-mstar+9
    allpost_super2.sfr0 = sfr0
    allpost_super2.Z = Z
    allpost_super2.tau = tau
    allpost_super2.chi2 = post_super2.chi2
    allpost_super2.scale = post_super2.scale

    for ii = 0, ngal-1 do begin
       temp = replicate({zobj: allised_super2[0].zobj, chi2: 0.0, chunkindx: 0L, $
         modelindx: 0L, ageindx: 0L, scale: 0.0},ndraw)
       temp.chunkindx = chunkindx[*,ii]
       temp.modelindx = modelindx[*,ii]
       temp.ageindx = ageindx[*,ii]
       temp.scale = allpost_super2[ii].scale
          
       postmodel1 = isedfit_restore(isedfit_paramfile,in_isedfit=temp,/fnu,$
         supergrid_paramfile=supergrid_paramfile,thissupergrid=2,$
         isedfit_dir=isedfit_dir,montegrids_dir=montegrids_dir)
       postmodel1.flux = postmodel1.flux*1D29 ; convert to uJy
       if ii eq 0 then postmodel_super2 = postmodel1 else postmodel_super2 = [[postmodel_super2],[postmodel1]]

;      plot, postmodel[0].wave/1D4, postmodel[0].flux, /xlog, xr=[1,7], yr=[-0.05,0.7], ps=3
;      for bb = 1, ndraw-1 do djs_oplot, postmodel[bb].wave/1D4, postmodel[bb].flux, ps=3
          
; get the posterior for beta; for broadband filters beta can be
; derived from the effective wavelengths of the filters and the color
; (see, e.g., Dunlop+12, Section 3):
;
;   beta = F160W-[ch1]/(-2.5*alog10(F160W_wave/[ch2]_wave)) - 2

       coeff = -2.5*alog10(betawave[0]/betawave[1])
       for jj = 0L, ndraw-1 do begin
          flam = postmodel1[jj].flux*1D-29*light/postmodel1[jj].wave^2
          mm = reform(k_project_filters(k_lambda_to_edges(postmodel1[jj].wave),$
            flam,filterlist=['clash_wfc3_f160w.par','spitzer_irac_ch1.par'],$
            /silent))
;         mmlam = reform(mm)*10^(-0.4*48.6)*light/betawave^2
          color = -2.5*alog10(mm[0]/mm[1])
          beta1 = color/coeff-2.0

;         djs_plot, postmodel1[0].wave/1D4, flam, xr=[1,6], xsty=1
;         djs_oplot, betawave/1D4, mmlam, psym=8, color='green'
          allpost_super2[ii].beta[jj] = beta1
       endfor
    endfor 
    
return
end    

pro z11_plots
; jm13jul30siena - make plots for the z11 galaxy Spitzer proposal

    common com_rebuild

    prefix = 'z11'
    rootpath = getenv('IM_PROJECTS_DIR')+'/clash/z11/'
    isedfit_dir = rootpath
    montegrids_dir = rootpath+'montegrids/'
    isedfit_paramfile = isedfit_dir+prefix+'_paramfile.par'
    sfhgrid_paramfile = isedfit_dir+prefix+'_sfhgrid.par'
    supergrid_paramfile = isedfit_dir+prefix+'_supergrid.par'

    filters = z11_filterlist(nice=nice_filters,pivotwave=pivotwave,width=width)
    nfilt = n_elements(filters)

; restore the posteriors    
    if n_elements(allpost_super2) eq 0 then rebuild_posteriors

    prob0 = exp(-0.5*(allpost_super2[0].chi2-allised_super2[1].chi2))
    prob1 = exp(-0.5*(allpost_super2[1].chi2-allised_super2[1].chi2))
    
    these0 = where(prob0 gt 0.75)
    these1 = where(prob1 gt 0.75)
    
; --------------------------------------------------
; posterior plots
    logmu = alog10(7.0) ; magnification for JD2
;   qquant = [0.5,0.25,0.75]
    qquant = [0.5,0.05,0.95]
;   qquant = [0.5,0.5-errorf(1.0/sqrt(2))/2, 0.5+errorf(1.0/sqrt(2))/2]
    
    psfile = rootpath+'z11_posteriors.ps'
    im_plotconfig, 5, pos, psfile=psfile, xspace=0.2, yspace=1.0, $
      xmargin=[0.8,0.3], height=2.3*[1,1], width=3*[1,1], charsize=1.6

; stellar mass
    quant0 = weighted_quantile(10^(allpost_super2[0].mstar[these0]-logmu),prob1[these0],quant=qquant)
    quant1 = weighted_quantile(10^(allpost_super2[1].mstar[these1]-logmu),prob1[these1],quant=qquant)
    best0 = where(10^(allpost_super2[0].mstar[these0]-logmu) gt quant0[1] and $
      10^(allpost_super2[0].mstar[these0]-logmu) lt quant0[2])
    best1 = where(10^(allpost_super2[1].mstar[these1]-logmu) gt quant1[1] and $
      10^(allpost_super2[1].mstar[these1]-logmu) lt quant1[2])

    splog, 'Stellar mass (5 hours) ', alog10(quant0[0]), alog10(quant0[1:2]/quant0[0])
    splog, 'Stellar mass (50 hours) ', alog10(quant1[0]), alog10(quant1[1:2]/quant1[0])
    
    bin = 0.05
    render_postplot, 10^(allpost_super2[1].mstar[these1]-logmu), pos[*,0], xrange=[1D7,1D10], $
      ytitle='', /logbins, binsize=bin, /xlog, $
      xtitle='Stellar Mass [(\mu/7)^{-1} M'+sunsymbol()+']', line=0, color_outline='black'
;   render_postplot, 10^(allpost_super2[1].mstar[these1[best1]]-logmu), pos[*,0], /overplot, $
;     /logbins, binsize=bin, color_fill='navy'
    render_postplot, 10^(allpost_super2[0].mstar[these0]-logmu), pos[*,0], /overplot, $
      binsize=bin, color_outline='black', /nofill, /logbins, line=1, thick=6

;   im_legend, 
    
; sSFR    
    quant0 = weighted_quantile(allpost_super2[0].sfrm[these0],prob0[these0],quant=qquant)
    quant1 = weighted_quantile(allpost_super2[1].sfrm[these1],prob1[these1],quant=qquant)

    best0 = where(allpost_super2[0].sfrm[these0] gt quant0[1] and $
      allpost_super2[0].sfrm[these0] lt quant0[2])
    best1 = where(allpost_super2[1].sfrm[these1] gt quant1[1] and $
      allpost_super2[1].sfrm[these1] lt quant1[2])

    print
    splog, 'sSFR (5 hours) ', 10^quant0[0], -(10^quant0[0]-10^quant0[1]), 10^quant0[2]-10^quant0[0]
    splog, 'sSFR (50 hours) ', 10^quant1[0], -(10^quant1[0]-10^quant1[1]), 10^quant1[2]-10^quant1[0]

    bin = 0.06
    render_postplot, allpost_super2[1].sfrm[these1], pos[*,1], /noerase, xrange=[-1.5,2.5], $
      xtitle='log_{10} sSFR [Gyr^{-1}]', bin=bin, line=0, color_outline='black'
    render_postplot, allpost_super2[0].sfrm[these0], pos[*,1], /overplot, bin=bin, $
      /nofill, color_outline='black', line=1, thick=6, /logbins, /xlog

;   bin = 0.05
;   render_postplot, 10^(allpost_super2[1].sfrm[these1]), pos[*,1], /noerase, xrange=[0.1,100], $
;     xtitle='sSFR [Gyr^{-1}]', bin=bin, line=0, color_outline='black', /logbins, /xlog
;   render_postplot, 10^allpost_super2[0].sfrm[these0], pos[*,1], /overplot, bin=bin, $
;     /nofill, color_outline='black', line=1, thick=6, /logbins, /xlog

; SFR-weighted age
    quant0 = weighted_quantile(allpost_super2[0].sfrage[these0]*1E3,prob0[these0],quant=qquant)
    quant1 = weighted_quantile(allpost_super2[1].sfrage[these1]*1E3,prob1[these1],quant=qquant)

    best0 = where(allpost_super2[0].sfrage[these0]*1E3 gt quant0[1] and $
      allpost_super2[0].sfrage[these0]*1E3 lt quant0[2])
    best1 = where(allpost_super2[1].sfrage[these1]*1E3 gt quant1[1] and $
      allpost_super2[1].sfrage[these1]*1E3 lt quant1[2])

    print
    splog, 'Age (5 hours) ', quant0[0], -(quant0[0]-quant0[1]), quant0[2]-quant0[0]
    splog, 'Age (50 hours) ', quant1[0], -(quant1[0]-quant1[1]), quant1[2]-quant1[0]
    
    bin = 15.0
    render_postplot, allpost_super2[1].sfrage[these1]*1E3, pos[*,2], xrange=[-20,450], /noerase, $
      xtitle='Age [SFR-weighted, Myr]', xtickinterval=100, bin=bin, /nonorm, line=0, color_outline='black'
    render_postplot, allpost_super2[0].sfrage[these0]*1E3, pos[*,2], bin=bin, /overplot, /nonorm, $
      /nofill, color_outline='black', line=1, thick=6

; beta
    quant0 = weighted_quantile(allpost_super2[0].beta[these0],prob0[these0],quant=qquant)
    quant1 = weighted_quantile(allpost_super2[1].beta[these1],prob1[these1],quant=qquant)

    best0 = where(allpost_super2[0].beta[these0] gt quant0[1] and $
      allpost_super2[0].beta[these0] lt quant0[2])
    best1 = where(allpost_super2[1].beta[these1] gt quant1[1] and $
      allpost_super2[1].beta[these1] lt quant1[2])

    print
    splog, 'beta (5 hours) ', quant0[0], -(quant0[0]-quant0[1]), quant0[2]-quant0[0]
    splog, 'beta (50 hours) ', quant1[0], -(quant1[0]-quant1[1]), quant1[2]-quant1[0]

    bin = 0.04
    render_postplot, allpost_super2[1].beta[these1], pos[*,3], xrange=[-3.0,-0.5], /noerase, $
      xtitle='UV Slope \beta', xtickinterval=0.5, bin=bin, /nonorm, line=0, color_outline='black'
    render_postplot, allpost_super2[0].beta[these0], pos[*,3], bin=bin, /overplot, /nonorm, $
      /nofill, color_outline='black', line=1, thick=6
    
; titles    
    xyouts, pos[0,0]-0.06, (pos[3,0]-pos[1,2])/2.0+pos[1,2], 'Marginalized Posterior Probability', $
      align=0.5, orientation=90, /norm, charsize=1.8
    
    im_plotconfig, psfile=psfile, /psclose, /pdf

; ---------------------------------------------------------------------------
; SED plot    
    psfile = rootpath+'z11_sed.ps'
    im_plotconfig, 0, pos, psfile=psfile, $
      height=5.0, margin=[1.2,0.5], width=7.0
;     height=4.5, margin=[1.0,1.5], width=6.0

    xtitle1 = textoidl('Wavelength \lambda (\mu'+'m)')
;   xtitle1 = textoidl('Observed-frame wavelength (\mu'+'m)')
    xtitle2 = textoidl('Rest-frame wavelength (\mu'+'m)')
    ytitle1 = textoidl('Flux F_{\nu} (\mu'+'Jy)')
    ytitle2 = textoidl('Magnitude (AB)')

;   scale = 1D30
;   yrange = [-0.5,6]
    yrange1 = [-0.07,0.8]
;   yrange = [29.5,23.0]
    xrange1 = [0.95,7.0]
;   xrange1 = [0.2,7.0]

    yrange2 = -2.5*alog10(yrange1>1D-9)+23.9
    
;   ticks1 = [0.3,1.2,4]
    ticks1 = [1.0,1.2,1.5,2,3,4,5,6]
;   ticks1 = [0.2,0.3,0.4,0.6,0.8,1,2,3,4,5,6]
    ticks2 = [0.02,0.04,0.06,0.1,0.2,0.3,0.4,0.6]

    yticks2 = [24.2,25,26,27,28.0]

    xrange2 = xrange1/(1.0+allised_super2[0].zobj)
    plot, [0], [0], /nodata, xrange=xrange1, yrange=yrange1, $
      xsty=1, ysty=1, xlog=1, position=pos, xtitle=xtitle1, $
      ytitle=ytitle1, xticks=n_elements(ticks1)-1, xtickv=ticks1, $
      xtickformat='(G0)', yminor=2
;   axis, /yaxis, ysty=1, ytitle='', yrange=yrange2, $
;     yticks=n_elements(yticks2)-1, ytickv=yticks2, ytickformat='(G0)'
;   xyouts, (pos[2]-pos[0])/2+pos[0], pos[3]+0.07, textoidl(ytitle2), $
;     align=0.5, /norm

    djs_oplot, 10^!x.crange, [0,0], line=5, thick=3
    
;   plot, [0], [0], /nodata, xrange=xrange1, yrange=yrange, $
;     xsty=5, ysty=5, xlog=1, position=pos
;   polyfill, [postmodel.wave/1D4,reverse(postmodel.wave/1D4)], $
;     [min(postmodel.flux,dim=2),reverse(max(postmodel.flux,dim=2))], $
;     /data, color=im_color('pale turquoise'), noclip=0, /fill
;   plot, [0], [0], /nodata, /noerase, xrange=xrange1, yrange=yrange, $
;     xsty=9, ysty=1, xtitle=xtitle1, ytitle=ytitle1, xlog=1, $
;     position=pos, xticks=n_elements(ticks1)-1, xtickv=ticks1, $
;     xtickformat='(G0)', yminor=2
;   axis, /xaxis, xsty=1, xtitle='', xrange=xrange2, xlog=1, $
;     xticks=n_elements(ticks2)-1, xtickv=ticks2, xtickformat='(G0)'
;   xyouts, (pos[2]-pos[0])/2+pos[0], pos[3]+0.07, textoidl(xtitle2), $
;     align=0.5, /norm

    im_legend, 'MACS0647-JD2 (z=10.8)', /left, /top, box=0, margin=0
    im_legend, ['Observed IRAC (5 hours)','Estimated IRAC (50 hours)'], $
      /left, /top, box=0, symsize=[2.2,2.5], $
      psym=[16,14], color=['grey50','navy'], symthick=10, margin=1.4
;   im_legend, ['Model Fluxes','Observed Fluxes (5 hours)','Estimated Fluxes (50 hours)'], $
;     /left, /top, box=0, symsize=[2.0,2.2,2.5], $
;     psym=[6,16,14], color=['blue','dark red','dark green'], symthick=10, margin=1.4

; overplot the 50% most likely models for each set of observations
;   polyfill, [postmodel[these0,0].wave/1D4,reverse(postmodel[these0,0].wave/1D4)], $
;     [min(postmodel[these0,0].flux,dim=2),reverse(max(postmodel[these0,0].flux,dim=2))], $
;     /data, color=im_color('grey90'), noclip=0, /fill
    polyfill, [postmodel_super2[these1,1].wave/1D4,reverse(postmodel_super2[these1,1].wave/1D4)], $
      [min(postmodel_super2[these1,1].flux,dim=2),reverse(max(postmodel_super2[these1,1].flux,dim=2))], $
      /data, color=im_color('powder blue'), noclip=0, /fill
;   djs_oplot, postmodel_super2[these1[0],1].wave/1D4, max(postmodel_super2[these1,1].flux,dim=2), line=0
;   djs_oplot, postmodel_super2[these1[0],1].wave/1D4, min(postmodel_super2[these1,1].flux,dim=2), line=0

    djs_oplot, postmodel_super2[these0[0],0].wave/1D4, max(postmodel_super2[these0,0].flux,dim=2), line=1
    djs_oplot, postmodel_super2[these0[0],0].wave/1D4, min(postmodel_super2[these0,0].flux,dim=2), line=1

; overplot the best-fit models    
;   oplot, allmodel_super2[0].wave/1D4, allmodel_super2[0].flux, line=0, thick=5;, color=im_color('brown4')
    oplot, allmodel_super2[1].wave/1D4, allmodel_super2[1].flux, line=0, thick=5, color=im_color('black')

; best-fit photometry    
;   djs_oplot, pivotwave/1D4, allised[0].bestmaggies*10^(0.4*23.9), $
;     psym=symcat(6,thick=10), symsize=3.0, color=im_color('grey')
    
; current HST observations    
    col = 'firebrick'
    these = where(pivotwave/1D4 lt 3)
    dowid = pivotwave[these]/1D4 gt 1.3
    errmaggies = 1.0/sqrt(allised_super2[0].ivarmaggies)
    oploterror, pivotwave[these]/1D4, allised_super2[0].maggies[these]*10^(0.4*23.9), $
      width[these]/1D4/2*dowid, errmaggies[these]*10^(0.4*23.9), psym=symcat(16), $
      symsize=2.2, color=im_color(col), errcolor=im_color(col), errthick=8

;   oploterror, pivotwave[these[other]]/1D4, allised[0].maggies[these[other]]*10^(0.4*23.9), $
;     errmaggies[these[other]]*10^(0.4*23.9), psym=symcat(16), $
;     symsize=2.2, color=im_color(col), errcolor=im_color(col), errthick=8
;   djs_oplot, pivotwave[upper]/1D4, errmaggies[upper]*10^(0.4*23.9), $
;     psym=symcat(11,thick=6), symsize=2.4, color=im_color('dark green')

; current [ch1], [ch2] fluxes
    col = 'grey50'
    irac = where(pivotwave/1D4 gt 3)
    errmaggies = 1.0/sqrt(allised_super2[0].ivarmaggies)
    oploterror, pivotwave[irac[0]]/1D4, allised_super2[0].maggies[irac[0]]*10^(0.4*23.9), $
      width[irac[0]]/1D4/2, errmaggies[irac[0]]*10^(0.4*23.9)*0, psym=symcat(16), $
      symsize=2.2, color=im_color(col), errcolor=im_color(col), errthick=8
    oploterror, pivotwave[irac[0]]/1D4, allised_super2[0].maggies[irac[0]]*10^(0.4*23.9), $
      errmaggies[irac[0]]*10^(0.4*23.9), psym=symcat(16), $
      symsize=2.2, color=im_color(col), errcolor=im_color(col), errthick=8, $
      /hibar

    oploterror, pivotwave[irac[1]]/1D4, allised_super2[0].maggies[irac[1]]*10^(0.4*23.9), $
      width[irac[1]]/1D4/2, errmaggies[irac[1]]*10^(0.4*23.9), psym=symcat(16), $
      symsize=2.2, color=im_color(col), errcolor=im_color(col), errthick=8
    
; proposed [ch1], [ch2] fluxes
    col = 'navy'
    errmaggies = 1.0/sqrt(allised_super2[1].ivarmaggies)
    oploterror, pivotwave[irac]/1D4, allised_super2[1].maggies[irac]*10^(0.4*23.9), $
      width[irac]/1D4/2, errmaggies[irac]*10^(0.4*23.9), psym=symcat(14), $
      symsize=3.0, color=im_color(col), errcolor=im_color(col), errthick=8
    
    im_plotconfig, psfile=psfile, /psclose, /pdf

; ---------------------------------------------------------------------------
; some QAplots    
    prob0 = exp(-0.5*(allpost_super2[0].chi2-allised_super2[1].chi2))
    prob1 = exp(-0.5*(allpost_super2[1].chi2-allised_super2[1].chi2))

    psfile = rootpath+'qa_z11_plots.ps'
    im_plotconfig, 0, pos, psfile=psfile, height=5.0, xmargin=[1.2,0.3]

; probability distributions    
    djs_plot, [0], [0], /nodata, xrange=[0,1.05], yrange=[0,600], $
      xsty=1, ysty=1, xtitle='Probability \propto exp(-\chi^2/2)', $
      ytitle='Number of Models', position=pos
    im_plothist, prob0, bin=0.01, /over, /fill, fcolor='firebrick'
    im_plothist, prob1, bin=0.01, /over, /fill, fcolor='dodger blue'
    im_legend, ['Existing (5 hours)','Proposed (50 hours)'], $
      /left, /top, box=0, textcolor=['firebrick','dodger blue']

; probability vs age
    djs_plot, [0], [0], /nodata, xrange=[0,450.0], yrange=[0,1.05], $
      xsty=1, ysty=1, ytitle='Probability \propto exp(-\chi^2/2)', $
      xtitle='Age (SFR-weighted, Myr)', position=pos
;   djs_oplot, !x.crange, [1,1], line=1
    djs_oplot, allpost_super2[0].sfrage*1E3, prob0, psym=symcat(16), $
      color=im_color('firebrick')
    djs_oplot, allpost_super2[1].sfrage*1E3, prob1, psym=symcat(15), $
      color=im_color('dodger blue')
    im_legend, ['Existing (5 hours)','Proposed (50 hours)'], $
      /right, /bottom, box=0, textcolor=['firebrick','dodger blue']
    
    im_plotconfig, psfile=psfile, /psclose, /pdf

stop    
    
    
return
end
    
