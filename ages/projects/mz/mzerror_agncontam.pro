pro mzerror_agncontam, ps=ps
; jm10aug05ucsd - model the systematic effects of residual AGN
; contamination 

    mzpath = ages_path(/projects)+'mz/'
    pspath = ages_path(/papers)+'mz/FIG_MZ/'
    if keyword_set(ps) then suffix = '.ps' else suffix = '.eps'

; read the relevant samples    
    hiianc = read_mz_sample(/mzhii_ancillary,/sdss)
    hiimass = read_mz_sample(/mzhii_mass,/sdss)
    hiiispec = read_mz_sample(/mzhii_ispec,/sdss)
    hiioh = read_mz_sample(/mzhii_log12oh,/sdss,/nodust)
    hiiinfo = mzlz_grab_info(hiioh,hiianc,hiimass,/flux,/t04,/nolimit,/errcut)
    hiibptd = hiiispec[hiiinfo.indx].bpt_d
    
    agnanc = read_mz_sample(/mzagn_ancillary,/sdss)
    agnmass = read_mz_sample(/mzagn_mass,/sdss)
    agnispec = read_mz_sample(/mzagn_ispec,/sdss)
    agnoh = read_mz_sample(/mzagn_log12oh,/sdss,/nodust)
    agninfo = mzlz_grab_info(agnoh,agnanc,agnmass,/flux,/t04,/nolimit,/errcut)
    agnbptd = agnispec[agninfo.indx].bpt_d

; fit the MZ/HII relation
    mzfit = fit_mz_closedbox(hiiinfo.mass,hiiinfo.oh,hiiinfo.weight,$
      oh_err=hiiinfo.oh_err,binsize=binsize,minmass=8.5,$
      maxmass=maxmass,mingal=mingal,fit_minmass=fit_minmass,$
      fit_maxmass=fit_maxmass,verbose=verbose)

; compute the distance from the median relation for the HII and AGN
; samples
    hiiindx = where((hiiinfo.mass gt min(mzfit.bin_mass)) and $
      (hiiinfo.mass lt max(mzfit.bin_mass)+binsize/2.0))
    hiidist = hiiinfo.oh[hiiindx]-mz_closedbox(hiiinfo.mass[hiiindx],mzfit.coeff)

    agnindx = where((agninfo.mass gt min(mzfit.bin_mass)) and $
      (agninfo.mass lt max(mzfit.bin_mass)+binsize/2.0))
    agndist = agninfo.oh[agnindx]-mz_closedbox(agninfo.mass[agnindx],mzfit.coeff)

    hiibin = im_medxbin(hiibptd[hiiindx],hiidist,0.05,minx=0.0,minpts=100)
    agnbin = im_medxbin(agnbptd[agnindx],agndist,0.05,minx=0.0,minpts=100)

; match the AGN and HII binned points and compute a "mixture" model of
; 80% star-forming and 20% AGN    
    match, hiibin.xbin, agnbin.xbin, m1, m2
    nmatch = n_elements(m1) 
    matchbin = im_empty_structure(hiibin,select=['xbin','medy'],$
      ncopies=nmatch)
    matchbin.xbin = hiibin[m1].xbin
    for ii = 0, nmatch-1 do matchbin[ii].medy = im_weighted_mean([hiibin[m1[ii]].medy,$
      agnbin[m1[ii]].medy],weights=[0.8,0.2])
;   struct_print, matchbin
    
; make some plots
    ohrange = [8.25,9.3]
    massrange = [8.4,11.6]
    residrange = [-1,0.5]
    drange = [-0.05,1.4]
;   levels = [0.01,0.05,0.25,0.5,0.75,0.95,0.99]
    levels = [0.5,0.75,0.9,0.975]
    maxis = range(8.8,11.3,75)
    
    psfile = pspath+'mzerror_agncontam'+suffix
    im_plotconfig, 1, pos, psfile=psfile, xspace=0.2, $
      xmargin=[1.1,1.2], width=[4.0,4.0], height=4.0
; MZ relation
    mzplot_scatterplot, /sdss, hiiinfo.mass, hiiinfo.oh, /nogrey, $
      position=pos[*,0], xsty=1, ysty=1, xrange=massrange, yrange=ohrange, $
      xtitle=mzplot_masstitle(), ytitle=mzplot_ohtitle(/t04,/fluxcor), levels=levels, $
      xtickinterval=1, /nooutlier, ccolor=im_color('grey40',102)
;   djs_oplot, mzfit.bin_mass, mzfit.bin_oh, psym=6, $
;     symsize=1.3, color='blue', thick=4
;   djs_oplot, mzfit.bin_mass, mzfit.bin_oh, $
;     line=0, thick=10, color=im_color('orange red',10)
    djs_oplot, maxis, mz_closedbox(maxis,mzfit.coeff), thick=7, $
      line=0;, color=im_color('orange red',101)
;   djs_oplot, mzfit.bin_mass, poly(mzfit.bin_mass,mzfit.coeff), $
;     line=0, thick=10, color=im_color('goldenrod',10)

    mzplot_scatterplot, /sdss, agninfo.mass, agninfo.oh, $
      position=pos[*,0], /noerase, xsty=5, ysty=5, $
      xrange=massrange, yrange=ohrange, /nooutlier, $
      /nogrey, ccolor=im_color('sienna',11), cthick=4, levels=levels
; money plot
    mzplot_scatterplot, /sdss, hiibptd[hiiindx], hiidist, /nogrey, $
      position=pos[*,1], /noerase, xsty=1, ysty=1, ytickinterval=0.5, $
      xrange=drange, yrange=residrange, xtickinterval=0.5, /nooutlier, $
      xtitle='D (dex)', ytitle='', ytickname=replicate(' ',10), levels=levels, $
      ccolor=im_color('grey40',102)
    axis, /yaxis, ysty=1, ytitle='', yrange=residrange, ytickinterval=0.5
    xyouts, pos[2,1]+0.08, (pos[3,1]-pos[1,1])/2.0+pos[1,1], $
      'Distance from M-Z (dex)', /norm, align=0.5, orientation=270

    mzplot_scatterplot, /sdss, agnbptd[agnindx], agndist, $
      position=pos[*,1], /noerase, xsty=5, ysty=5, $
      xrange=drange, yrange=residrange, /nooutlier, $
      /nogrey, ccolor=im_color('sienna',11), cthick=4, levels=levels

    djs_oplot, hiibin.xbin, hiibin.medy, line=0, $
      color=im_color('orange red',10), thick=8
    djs_oplot, agnbin.xbin, agnbin.medy, line=5, $
      color='navy', thick=10
    djs_oplot, matchbin.xbin, matchbin.medy, line=3, $
      color='red', thick=10
;   djs_oplot, hiiispec.bpt_d, hiidist, psym=3, color='red'
    im_plotconfig, /psclose

stop
return
end
    
