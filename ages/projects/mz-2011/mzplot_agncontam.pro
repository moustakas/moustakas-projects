pro mzplot_agncontam, ps=ps
; jm10aug05ucsd - model the systematic effects of residual AGN
; contamination

    mzpath = mz_path()
    qapath = mzpath+'qaplots/'
    if keyword_set(ps) then begin
       pspath = qapath
       suffix = '.ps'
    endif else begin
       pspath = mz_path(/paper)
       suffix = '.eps'
    endelse

; read the star-forming sample
    hiianc = read_mz_sample(/mzhii_ancillary,/sdss)
    hiimass = read_mz_sample(/mzhii_mass,/sdss)
    hiiispec = read_mz_sample(/mzhii_ispec,/sdss)
    hiioh = read_mz_sample(/mzhii_log12oh,/sdss,/nodust)
    hiiinfo = mzlz_grab_info(hiioh,hiianc,hiimass,/flux,/t04,/nolimit,/errcut)
    hiibptd = hiiispec[hiiinfo.indx].bpt_d

; for the full AGN sample...
    allagnanc = read_mz_sample(/mzagn_ancillary,/sdss)
    allagnmass = read_mz_sample(/mzagn_mass,/sdss)
    allagnispec = read_mz_sample(/mzagn_ispec,/sdss)
    allagnoh = read_mz_sample(/mzagn_log12oh,/sdss,/nodust)

; ...but divide it into strong AGN and potential composite AGN, as
; classified by Yan
    iscomp = where(allagnispec.yan_agn eq 0 and allagnispec.bpt_agn eq 1,ncomp)
    isagn = where(allagnispec.yan_agn eq 1)

    agnanc = allagnanc[isagn]
    agnmass = allagnmass[isagn]
    agnispec = allagnispec[isagn]
    agnoh = allagnoh[isagn]

    companc = allagnanc[iscomp]
    compmass = allagnmass[iscomp]
    compispec = allagnispec[iscomp]
    compoh = allagnoh[iscomp]
    
    agninfo = mzlz_grab_info(agnoh,agnanc,agnmass,/flux,/t04,/nolimit,/errcut)
    agnbptd = agnispec[agninfo.indx].bpt_d

    compinfo = mzlz_grab_info(compoh,companc,compmass,/flux,/t04,/nolimit,/errcut)
    compbptd = compispec[compinfo.indx].bpt_d
    
;; ...but divide it into Kewley AGN and SF/AGN
;    class = iclassification(agnispec[agninfo.indx],$
;      snrcut_class=0.0,ratios=agnratios,doplot=1)
;    isagn = where(strtrim(agnratios.final_class eq 'AGN'),nisagn)
;    issfagn = where(strtrim(agnratios.final_class eq 'SF/AGN'),nissfagn)

; read the previously fitted MZ relation
    mzfit = mrdfits(mzpath+'mzlocal_sdss_fluxcor_t04.fits.gz',1)
;   mzfit = fit_mz_closedbox(hiiinfo.mass,hiiinfo.oh,hiiinfo.weight,$
;     oh_err=hiiinfo.oh_err,binsize=binsize,minmass=8.5,$
;     maxmass=maxmass,mingal=mingal,fit_minmass=fit_minmass,$
;     fit_maxmass=fit_maxmass,verbose=verbose)

    binsize = 0.1
    
; compute the distance from the MZ relation for the various subsamples 
    hiiindx = where((hiiinfo.mass gt min(mzfit.bin_mass)) and $
      (hiiinfo.mass lt max(mzfit.bin_mass)+binsize/2.0))
    hiidist = hiiinfo.oh[hiiindx]-mz_closedbox(hiiinfo.mass[hiiindx],mzfit.coeff)

    agnindx = where((agninfo.mass gt min(mzfit.bin_mass)) and $
      (agninfo.mass lt max(mzfit.bin_mass)+binsize/2.0))
    agndist = agninfo.oh[agnindx]-mz_closedbox(agninfo.mass[agnindx],mzfit.coeff)

    compindx = where((compinfo.mass gt min(mzfit.bin_mass)) and $
      (compinfo.mass lt max(mzfit.bin_mass)+binsize/2.0))
    compdist = compinfo.oh[compindx]-mz_closedbox(compinfo.mass[compindx],mzfit.coeff)

    hiibin = im_medxbin(hiibptd[hiiindx],hiidist,0.05,minx=0.0,minpts=100)
    agnbin = im_medxbin(agnbptd[agnindx],agndist,0.05,minx=0.0,minpts=100)
    compbin = im_medxbin(compbptd[compindx],compdist,0.05,minx=0.0,minpts=100)

; match the AGN and HII binned points and compute a "mixture" model of
; 80% star-forming and 20% AGN    
    match, hiibin.xbin, compbin.xbin, m1, m2
    niceprint, hiibin[m1].xbin, compbin[m2].xbin, hiibin[m1].medy, compbin[m2].medy, $
      hiibin[m1].medy-compbin[m2].medy, 0.5*(hiibin[m1].medy-compbin[m2].medy)
    
    nmatch = n_elements(m1) 
    matchbin = im_empty_structure(hiibin,select=['xbin','medy'],$
      ncopies=nmatch)
    matchbin.xbin = hiibin[m1].xbin
    for ii = 0, nmatch-1 do matchbin[ii].medy = im_weighted_mean([hiibin[m1[ii]].medy,$
      compbin[m2[ii]].medy],weights=[0.8,0.2])
;   struct_print, matchbin
    
; make some plots
    ohrange = [8.2,9.25]
    massrange = [9.0,11.5]
    residrange = [-1,0.5]
    drange = [-0.05,1.4]
;   levels = [0.01,0.05,0.25,0.5,0.75,0.95,0.99]
    levels = [0.25,0.5,0.75,0.9,0.975]
    maxis = range(9.2,11.3,75)

    npix1 = 29
    npix2 = 21
    
    psfile = pspath+'mzerror_agncontam'+suffix
    im_plotconfig, 1, pos, psfile=psfile, xspace=0.2, $
      xmargin=[1.1,1.2], width=[4.0,4.0], height=4.0
; MZ relation
    mzplot_scatterplot, /sdss, hiiinfo.mass, hiiinfo.oh, /nogrey, $
      position=pos[*,0], xsty=1, ysty=1, xrange=massrange, yrange=ohrange, $
      xtitle=mzplot_masstitle(), ytitle=mzplot_ohtitle(/t04,/fluxcor), $
      levels=levels, xnpix=npix1, ynpix=npix1, $
      xtickinterval=1, ccolor=djs_icolor('grey'), cthick=3, /nooutlier
;   djs_oplot, mzfit.bin_mass, mzfit.bin_oh, psym=6, $
;     symsize=1.3, color='blue', thick=4
;   djs_oplot, mzfit.bin_mass, mzfit.bin_oh, $
;     line=0, thick=10, color=im_color('orange red red',10)
;   djs_oplot, maxis, mz_closedbox(maxis,mzfit.coeff), thick=7, line=0
;   djs_oplot, mzfit.bin_mass, poly(mzfit.bin_mass,mzfit.coeff), $
;     line=0, thick=10, color=im_color('goldenrod',10)

    mzplot_scatterplot, /sdss, agninfo.mass, agninfo.oh, $
      position=pos[*,0], /noerase, xsty=5, ysty=5, $
      xrange=massrange, yrange=ohrange, nooutlier=1, $
      /nogrey, ccolor=im_color('orange'), cthick=3, $
      levels=levels, xnpix=npix2, ynpix=npix2, cline=0, $
      outcolor=im_color('orange')

    mzplot_scatterplot, /sdss, compinfo.mass, compinfo.oh, $
      position=pos[*,0], /noerase, xsty=5, ysty=5, $
      xrange=massrange, yrange=ohrange, nooutlier=1, $
      /nogrey, ccolor=im_color('blue'), cthick=3, $
      levels=levels, xnpix=npix1, ynpix=npix1, cline=0, $
      outcolor=im_color('blue')

; money plot
    mzplot_scatterplot, /sdss, hiibptd[hiiindx], hiidist, /nogrey, $
      position=pos[*,1], /noerase, xsty=1, ysty=1, ytickinterval=0.5, $
      xrange=drange, yrange=residrange, xtickinterval=0.5, $
      xtitle='D (dex)', ytitle='', ytickname=replicate(' ',10), levels=levels, $
      ccolor=djs_icolor('grey'), xnpix=npix1, ynpix=npix1, cthick=3, /nooutlier
    axis, /yaxis, ysty=1, ytitle='', yrange=residrange, ytickinterval=0.5
    xyouts, pos[2,1]+0.08, (pos[3,1]-pos[1,1])/2.0+pos[1,1], $
      'Distance from M-Z (dex)', /norm, align=0.5, orientation=270

    mzplot_scatterplot, /sdss, agnbptd[agnindx], agndist, $
      position=pos[*,1], /noerase, xsty=5, ysty=5, $
      xrange=drange, yrange=residrange, nooutlier=1, $
      /nogrey, ccolor=im_color('orange'), cthick=3, $
      levels=levels, xnpix=npix2, ynpix=npix2, cline=0, $
      outcolor=im_color('orange')

    mzplot_scatterplot, /sdss, compbptd[compindx], compdist, $
      position=pos[*,1], /noerase, xsty=5, ysty=5, $
      xrange=drange, yrange=residrange, nooutlier=1, $
      /nogrey, ccolor=im_color('blue'), cthick=3, $
      levels=levels, xnpix=npix1, ynpix=npix1, cline=0, $
      outcolor=im_color('blue')

    djs_oplot, hiibin.xbin, hiibin.medy, line=0, color=im_color('black'), thick=8
    djs_oplot, compbin.xbin, compbin.medy, line=5, color=im_color('dodger blue'), thick=8
    djs_oplot, agnbin.xbin, agnbin.medy, line=3, color=im_color('firebrick'), thick=8

;   djs_oplot, matchbin.xbin, matchbin.medy, line=3, color='red', thick=10
;   djs_oplot, hiiispec.bpt_d, hiidist, psym=3, color='red'
    im_plotconfig, /psclose

stop
return
end
    
