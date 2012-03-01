pro render_postplot, xx, pos, xrange=xrange, yrange=yrange, $
  xtitle=xtitle, ytitle=ytitle, binsize=binsize, noerase=noerase, $
  xtickinterval=xtickinterval, nomedian=nomedian, monte=monte

    yrange = [0,1.1]
;   yrange = [0,1.05]

    if (n_elements(xrange) eq 0) then xrange = minmax(xx)*[0.9,1.1]
    if (n_elements(binsize) eq 0) then $
      binsize = (xrange[1]-xrange[0])/ceil(0.3*sqrt(n_elements(xx)))

    plot, [0], [0], xsty=5, ysty=5, /nodata, position=pos, $
      yrange=yrange, xrange=xrange, noerase=noerase
    im_plothist, xx, bin=binsize, /peak, /overplot, /fill, $
      fcolor=im_color('grey80'), xhist, yhist
;   if keyword_set(nomedian) eq 0 then $
;     djs_oplot, median(xx)*[1,1], !y.crange, line=5, thick=8
    plot, [0], [0], xsty=1, ysty=1, /noerase, /nodata, position=pos, $
      yrange=yrange, xrange=xrange, ytitle=ytitle, $
      xtitle=textoidl(xtitle), xtickinterval=xtickinterval, $
      ytickname=replicate(' ',10)
    if n_elements(monte) ne 0 then begin
       im_plothist, monte, bin=binsize*2, mxhist, myhist, /noplot
       im_plothist, monte, bin=binsize*2, /overplot, line=1, $
         normfactor=max(myhist)/(yrange[1]*0.95)
;        normfactor=max(myhist)/max(yhist)/1.2
    endif

return
end

pro santorini_rebuild_posteriors

    common com_santorini, allpost, allised, postmodel
    
    datapath = clash_path(/santorini)
    isedpath = datapath+'isedfit/'

    isedfit_sfhgrid_dir = datapath+'montegrids/'
    sfhgrid_paramfile = getenv('CLASH_DIR')+'/santorini/santorini_sfhgrid.par'
    params = yanny_readone(sfhgrid_paramfile)
    nmodel = params.nmonte
    ndraw = isedfit_ndraw()

    cat = read_santorini()

    supergrid = [1,2,3,4]
    super = get_santorini_supergrid(supergrid,nsuper=nsuper)
;   super = super[0:2] & nsuper = 3
    struct_print, super
    
    allpost = {mstar: fltarr(ndraw), sfrage: fltarr(ndraw), Z: fltarr(ndraw), $
      age: fltarr(ndraw), tau: fltarr(ndraw), sfr0: fltarr(ndraw), $
      av: fltarr(ndraw)};, bigsfr: fltarr(nmodel), bigsfrage: fltarr(nmodel), $
;     bigmass: fltarr(nmodel)}
    allpost = replicate(allpost,nsuper)

    for ii = 0, nsuper-1 do begin
       paramfile = isedpath+'santorini_supergrid0'+strtrim(ii+1,2)+'_isedfit.par'
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

;      allpost[ii].bigsfr = alog10(bigsfr)
;      allpost[ii].bigmass = alog10(bigmass)
;      allpost[ii].bigsfrage = bigsfrage

; reconstruct the posterior models, but just for the fiducial grid 
       if ii eq 1 then begin
          temp = replicate({zobj: cat.z, chi2: 0.0, chunkindx: 0L, $
            modelindx: 0L, ageindx: 0L, scale: 0.0},ndraw)
          temp.chunkindx = chunkindx
          temp.modelindx = modelindx
          temp.ageindx = ageindx
          temp.scale = post.scale
;         temp.chi2 = post.chi2

          splog, paramfile
          postmodel = isedfit_restore(paramfile,in_isedfit=temp,$
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

    filters = clash_filterlist(nice=nice_filters,/useirac,weff=weff)
    ndraw = isedfit_ndraw()

; gather the photometry and restore the results
    cat = read_santorini()
    nobj = n_elements(cat)

    mu = 17.0
    logmu = alog10(mu)
    
    if n_elements(allpost) eq 0 then santorini_rebuild_posteriors

    niceprint, allised.mass_50-logmu, allised.sfr_50-logmu, allised.sfrage_50*1E3, $
      allised.mass-logmu, allised.sfr-logmu, allised.sfrage*1E3
    
    splog, weighted_quantile(allpost[1].sfrage*1E3,quant=0.95)
    splog, getredshift(getage(9.6)-0.220)

; ---------------------------------------------------------------------------
; SED plot for the paper
    paramfile = isedpath+'santorini_supergrid02_isedfit.par'
    model = isedfit_restore(paramfile,ised,iopath=isedpath,$
      isedfit_sfhgrid_dir=isedfit_sfhgrid_dir)

    psfile = datapath+'santorini_supergrid02_sed.eps'
    im_plotconfig, 0, pos, psfile=psfile, height=5.0, ymargin=[1.0,1.0]

    xtitle1 = textoidl('Observed-Frame Wavelength (\AA)')
    xtitle2 = textoidl('Rest-Frame Wavelength (\AA)')
    ytitle1 = textoidl('Apparent AB Magnitude')

    yrange = [29.0,24.0]
    xrange1 = [1900,70000]

    ticks1 = [3000,12000,40000]
    ticks2 = [250,500,1000,2000,4000]

    xrange2 = xrange1/(1.0+ised.zobj)
    plot, [0], [0], /nodata, xrange=xrange1, yrange=yrange, $
      xsty=9, ysty=1, xtitle=xtitle1, ytitle=ytitle1, xlog=1, $
      position=pos, xtickformat='(I0)', xticks=n_elements(ticks1)-1, xtickv=ticks1
    axis, /xaxis, xsty=1, xtitle=xtitle2, xrange=xrange2, xlog=1, $
      xticks=n_elements(ticks2)-1, xtickv=ticks2

    nww = 100
    ww = shuffle_indx(ndraw,num=nww)

;    for ii = 0, nww-1 do begin
;       djs_oplot, postmodel[ww[ii]].wave, postmodel[ww[ii]].flux, $
;         line=0, thick=1, color=im_color('grey80')
;;      cc = get_kbrd(1)
;    endfor

;   polyfill, [postmodel[0].wave,reverse(postmodel[0].wave)], $
;     [min(postmodel.flux,dim=2),reverse(max(postmodel.flux,dim=2))], $
;     /data, color=im_color('tan'), noclip=0, /fill
;   im_legend, [cat[ii].galaxy,'z_{phot} = '+string(cat[0].z,format='(F4.2)')], $
;     /left, /top, box=0, margin=0, charsize=1.7
    oplot, model[0].wave, model[0].flux, line=0
    
    mab = maggies2mag(ised[0].maggies,ivar=ised[0].ivarmaggies,$
      lomagerr=loerr,himagerr=hierr,magnsigma=mabupper,nsigma=2.0)
    good = where(mab gt -99.0,ngood,comp=upper,ncomp=nupper)

    djs_oplot, weff[good], -2.5*alog10(ised[0].bestmaggies[good]), $
      psym=symcat(6,thick=6), symsize=2.5, color=im_color('steel blue')

    if (ngood ne 0L) then begin
       oploterror, weff[good], mab[good], loerr[good], psym=symcat(16), $
         symsize=2.0, color=im_color('firebrick'), $
         errcolor=im_color('firebrick'), errthick=8, /lobar
       oploterror, weff[good], mab[good], hierr[good], psym=3, $
         symsize=2.0, color=im_color('firebrick'), $
         errcolor=im_color('firebrick'), errthick=8, /hibar
    endif

    if (nupper ne 0) then begin
       djs_oplot, weff[upper], mabupper[upper], psym=symcat(11,thick=6), $
         symsize=3.0, color=im_color('grey20')
    endif

    im_plotconfig, /psclose, psfile=psfile, /pdf, /pskeep

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
; posterior distributions plots

; supergrid01    
    monte = mrdfits(isedfit_sfhgrid_dir+'sfhgrid01/fsps/calzetti/chab_montegrid.fits.gz',1)

    psfile = datapath+'santorini_supergrid01_posteriors.eps'
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

; supergrid02
    monte = mrdfits(isedfit_sfhgrid_dir+'sfhgrid02/fsps/chab_montegrid.fits.gz',1)

    psfile = datapath+'santorini_supergrid02_posteriors.eps'
    im_plotconfig, 5, pos, psfile=psfile, xspace=0.3, yspace=1.0, $
      xmargin=[1.0,0.5], height=3.5*[1,1], width=4*[1,1], charsize=1.8

    render_postplot, allpost[1].mstar-logmu, pos[*,0], xrange=[6.5,9.5], $
      ytitle='Probability', xtitle='log (M_{*}/M'+sunsymbol()+')', $
      xtickinterval=1
;     monte=allpost[1].bigmass+median(allpost[1].mstar-logmu)+1.5

    render_postplot, allpost[1].sfr0-logmu, pos[*,1], /noerase, xrange=[-0.5,0.3], $ ; xrange=minmax(allpost[1].sfr0-logmu)*[0.99,1.01], $
      xtitle='log (SFR /M'+sunsymbol()+' yr^{-1})';, monte=allpost[1].bigsfr+10
    
    render_postplot, allpost[1].tau, pos[*,2], /noerase, xrange=[-0.2,1.2], $
      ytitle='Probability', xtitle='\tau (Gyr)', xtickinterval=0.5;monte=monte.tau, $
    
    render_postplot, allpost[1].sfrage*1E3, pos[*,3], xrange=[-10,500], /noerase, $
      xtitle='<t>_{SFR} (Myr)', xtickinterval=100
;   im_plothist, allpost[1].sfrage*1E3, /cumu, /overplot, color='red'

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

return
end
