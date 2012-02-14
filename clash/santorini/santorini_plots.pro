;pro rebuild_posteriors
;
;    common com_santorini, allpost, postmodel
;    
;    ngrid = 3
;    ndraw = isedfit_ndraw()
;    allpost = {mstar: fltarr(ndraw), sfrage: fltarr(ndraw), 
;
;    for
;    
;    paramfile = isedpath+'santorini_supergrid01_isedfit.par'
;
;    
;    
;
;return
;end    


pro santorini_plots, supergrid, models=models, isedfit=isedfit, $
  qaplot=qaplot, clobber=clobber
; jm11nov08ucsd - build plots for the paper

    datapath = clash_path(/santorini)
    isedpath = datapath+'isedfit/'

    isedfit_sfhgrid_dir = datapath+'montegrids/'
    sfhgrid_paramfile = getenv('CLASH_DIR')+'/santorini/santorini_sfhgrid.par'

    filters = clash_filterlist(nice=nice_filters,/useirac,weff=weff)
    ndraw = isedfit_ndraw()

; gather the photometry and restore the results
    cat = read_santorini()
    nobj = n_elements(cat)

    common com_post, mstar, post, age, Z, tau, sfr0, b100, $
      av, sfrage, chunkindx, modelindx, ageindx, postmodel
    paramfile = isedpath+'santorini_supergrid01_isedfit.par'

    if (n_elements(mstar) eq 0L) then begin
       mstar = isedfit_reconstruct_posterior(paramfile,post=post,$
         isedfit_sfhgrid_dir=isedfit_sfhgrid_dir,iopath=isedpath,$
         age=age,Z=Z,tau=tau,sfr0=sfr0,b100=b100,av=av,sfrage=sfrage,$
         chunkindx=chunkindx,modelindx=modelindx,indxage=ageindx)

; reconstruct the posterior models       
       temp = replicate({zobj: cat.z, chi2: 0.0, chunkindx: 0L, $
         modelindx: 0L, ageindx: 0L, scale: 0.0},ndraw)
       temp.chunkindx = chunkindx
       temp.modelindx = modelindx
       temp.ageindx = ageindx
       temp.scale = post.scale
       postmodel = isedfit_restore(paramfile,in_isedfit=temp,$
         isedfit_sfhgrid_dir=isedfit_sfhgrid_dir,iopath=isedpath)
;      plot, postmodel[0].wave, postmodel[0].flux, /xlog, xr=[1E4,4E4], yr=[30,20], ps=3
;      for ii = 1, ndraw-1 do djs_oplot, postmodel[ii].wave, postmodel[ii].flux, ps=3
    endif

; ---------------------------------------------------------------------------
; SED plot for the paper
    paramfile = isedpath+'santorini_supergrid01_isedfit.par'
    model = isedfit_restore(paramfile,ised,iopath=isedpath,$
      isedfit_sfhgrid_dir=isedfit_sfhgrid_dir)

    psfile = datapath+'santorini_sed.eps'
    im_plotconfig, 0, pos, psfile=psfile, height=5.0, ymargin=[1.0,1.0]
;   im_plotfaves

    xtitle1 = textoidl('Observed-Frame Wavelength (\AA)')
    xtitle2 = textoidl('Rest-Frame Wavelength (\AA)')
    ytitle1 = textoidl('Apparent AB Magnitude')

    yrange = [30.0,24.0]
    xrange1 = [1900,70000]

    ticks1 = [4000,12000,30000]
    ticks2 = [250,500,1000,2000,4000]

    xrange2 = xrange1/(1.0+ised[0].zobj)
    plot, [0], [0], /nodata, xrange=xrange1, yrange=yrange, $
      xsty=9, ysty=1, xtitle=xtitle1, ytitle=ytitle1, xlog=1, $
      position=pos, xtickformat='(I0)', xticks=n_elements(ticks1)-1, xtickv=ticks1
    axis, /xaxis, xsty=1, xtitle=xtitle2, xrange=xrange2, xlog=1, $
      xticks=n_elements(ticks2)-1, xtickv=ticks2
;   for ii = 0, 75 do begin
;      djs_oplot, postmodel[ii].wave, postmodel[ii].flux, $
;        line=0, thick=1, color=im_color('grey80')
;      cc = get_kbrd(1)
;   endfor
;   prob = exp(-0.5*(post.chi2-ised.chi2))
;   prob = prob/total(prob,/double)
;   ww = where(prob gt weighted_quantile(prob,quant=0.25
    nww = 100
;   ww = shuffle_indx(ndraw,num=nww)
    ww = [1017,4,545,153,854,1841,438,1868,415,605,1438,1141,595,1558,1926,1241,1206,108,163,1597,1454,$
      199,313,1680,1108,1853,641,1764,1742,1739,1789,1727,1712,1872,1902,1313,57,1549,1270,556,1699,1039,$
      1251,112,191,1297,629,1961,1376,921,195,1968,859,1542,242,1927,621,1265,961,550,1358,282,1192,1380,$
      54,1803,979,337,1954,1867,483,1592,1488,38,235,27,1532,1303,1741,1304,568,35,1552,1472,797,372,239,$
      81,1426,880,1576,948,1814,1564,279,1170,179,1599,66,1041]
    for ii = 0, nww-1 do begin
       djs_oplot, postmodel[ww[ii]].wave, postmodel[ww[ii]].flux, $
         line=0, thick=1, color=im_color('grey80')
;      cc = get_kbrd(1)
    endfor
;   polyfill, [postmodel[0].wave,reverse(postmodel[0].wave)], $
;     [min(postmodel.flux,dim=2),reverse(max(postmodel.flux,dim=2))], $
;     /data, color=im_color('tan'), noclip=0, /fill
;   im_legend, [cat[ii].galaxy,'z_{phot} = '+string(cat[0].z,format='(F4.2)')], $
;     /left, /top, box=0, margin=0, charsize=1.7
    oplot, model[0].wave, model[0].flux, line=0
    
    used = where((ised[0].maggies gt 0.0) and $ ; used in the fitting
      (ised[0].ivarmaggies gt 0.0),nused)
    upper = where((ised[0].maggies le 0.0) and $ ; upper limit
      (ised[0].ivarmaggies gt 0.0),nupper)

    djs_oplot, weff[used], -2.5*alog10(ised[0].bestmaggies[used]), $
      psym=symcat(6,thick=6), symsize=2.5, color=im_color('steel blue')

    if (nused ne 0L) then begin
       mab = maggies2mag(ised[0].maggies[used],$
         ivar=ised[0].ivarmaggies[used],magerr=mab_err)
       oploterror, weff[used], mab, mab_err, psym=symcat(16), $
         symsize=2.0, color=im_color('firebrick'), $
         errcolor=im_color('firebrick'), errthick=8
    endif

    if (nupper ne 0) then begin
       mab = maggies2mag(1.0/sqrt(ised[0].ivarmaggies[upper]))
       djs_oplot, weff[upper], mab, psym=symcat(11,thick=6), $
         symsize=3.0, color=im_color('dark orchid4')
    endif

    im_plotconfig, /psclose, psfile=psfile, /pdf, /pskeep

stop    
    
; ---------------------------------------------------------------------------
; compare the posterior distributions on age 
    paramfile1 = isedpath+'santorini_supergrid01_isedfit.par'
    mstar1 = isedfit_reconstruct_posterior(paramfile1,sfrage=sfrage1,$
      isedfit_sfhgrid_dir=isedfit_sfhgrid_dir,iopath=isedpath)

    paramfile2 = isedpath+'santorini_supergrid02_isedfit.par'
    mstar2 = isedfit_reconstruct_posterior(paramfile2,sfrage=sfrage2,$
      isedfit_sfhgrid_dir=isedfit_sfhgrid_dir,iopath=isedpath)

    paramfile3 = isedpath+'santorini_supergrid03_isedfit.par'
    mstar3 = isedfit_reconstruct_posterior(paramfile3,sfrage=sfrage3,$
      isedfit_sfhgrid_dir=isedfit_sfhgrid_dir,iopath=isedpath)
        
    
    psfile = datapath+'santorini_sfrage.eps'
    im_plotconfig, 0, pos, psfile=psfile, height=5.0

    xr = [-10,500]
    yr = [0,1.1]
    bin = 50
    plot, [0], [0], xsty=5, ysty=5, /nodata, yrange=yr, $
      xrange=xr, position=pos, ytitle='', xtitle='', charsize=1.4
    im_plothist, sfrage1*1E3, bin=bin, /peak, /overplot, /fill, color=im_color('green');, fcolor=im_color('grey80')
    im_plothist, sfrage2*1E3, bin=bin, /peak, /overplot, /fill, color=im_color('orange')
    im_plothist, sfrage3*1E3, bin=bin, /peak, /overplot, /fill, color=im_color('red')
    plot, [0], [0], xsty=1, ysty=1, /noerase, /nodata, yrange=yr, $
      xrange=xr, position=pos, ytitle='', xtitle=textoidl('SFR-Weighted Age (Myr)'), $
      ytickname=replicate(' ',10), charsize=1.4, xtickinterval=0.5
    
    im_plotconfig, /psclose, psfile=psfile, /pdf, /pskeep

stop    
    

; ---------------------------------------------------------------------------
; plot the NxN distributions of parameters as an upper triangle
    this = 0 ; z=9.56 solution
    
    psfile = datapath+'santorini_manyd.eps'
    splog, 'Writing '+psfile
    manyd = transpose([[mstar[*,this]],[sfrage[*,this]],[Z[*,this]/0.019],[av[*,this]],[sfr0[*,this]]])
    label = textoidl(['log (M/M'+sunsymbol()+')','t_{w} (Gyr)','Z/Z'+sunsymbol(),'A_{V} (mag)','log (\psi/M'+sunsymbol()+' yr^{-1})'])
    
    im_manyd_scatterplot, fltarr(ndraw)+1, manyd, psfile, label=label, $
      axis_char_scale=1.4, /internal, outliers=1, $
      /nogrey, levels=errorf((dindgen(2)+1)/sqrt(2)), /upper
    spawn, 'ps2pdf '+psfile, /sh

stop    

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


return
end
