pro mzerror_agncontam, ps=ps
; jm10aug05ucsd - model the systematic effects of residual AGN
; contamination 

    mzpath = ages_path(/projects)+'mz/'
    pspath = ages_path(/papers)+'mz/FIG_MZ/'
    if keyword_set(ps) then suffix = '.ps' else suffix = '.eps'

; read the relevant samples    
;   hiianc = read_mz_sample(/mzhii_ancillary,/sdss)
    hiimass = read_mz_sample(/mzhii_mass,/sdss)
    hiiispec = read_mz_sample(/mzhii_ispec,/sdss)
    hiioh = read_mz_sample(/mzhii_log12oh,/sdss,/nodust)
    hiigood = where((strtrim(hiioh.r23branch,2) eq 'U') and $
      (hiiispec.bpt_d gt -900))
;   hiianc = hiianc[hiigood]
    hiimass = hiimass[hiigood]
    hiiispec = hiiispec[hiigood]
    hiioh = hiioh[hiigood]

    agnmass = read_mz_sample(/mzagn_mass,/sdss)
    agnispec = read_mz_sample(/mzagn_ispec,/sdss)
    agnoh = read_mz_sample(/mzagn_log12oh,/sdss,/nodust)
    agngood = where((strtrim(agnoh.r23branch,2) eq 'U') and $
      (agnispec.bpt_d gt -900))
    agnmass = agnmass[agngood]
    agnispec = agnispec[agngood]
    agnoh = agnoh[agngood]

; fit the MZ/HII relation
    mzfit = fit_mz(hiimass.mass_avg,hiioh.log12oh_t04,$
      minmass=8.5,binsize=binsize)

; compute the distance from the median relation for the HII and AGN
; samples
;   hiiindx = where((hiimass.mass_avg gt min(mzfit.bin_mass)))
    hiiindx = where((hiimass.mass_avg gt min(mzfit.bin_mass)) and $
      (hiimass.mass_avg lt max(mzfit.bin_mass)+binsize/2.0))
    hiidist = hiioh[hiiindx].log12oh_t04-poly(hiimass[hiiindx].mass_avg,mzfit.coeff)

;   agnindx = where((hiimass.mass_avg gt min(mzfit.bin_mass)))
    agnindx = where((agnmass.mass_avg gt min(mzfit.bin_mass)) and $
      (agnmass.mass_avg lt max(mzfit.bin_mass)+binsize/2.0))
    agndist = agnoh[agnindx].log12oh_t04-poly(agnmass[agnindx].mass_avg,mzfit.coeff)

    hiibin = im_medxbin(hiiispec[hiiindx].bpt_d,hiidist,0.05,minx=0.0,minpts=100)
    agnbin = im_medxbin(agnispec[agnindx].bpt_d,agndist,0.05,minx=0.0,minpts=100)

; match the AGN and HII binned points and compute a "mixture" model of
; 80% star-forming and 20% AGN    
    match, hiibin.xbin, agnbin.xbin, m1, m2
    nmatch = n_elements(m1) 
    matchbin = im_empty_structure(hiibin,select=['xbin','medy'],$
      ncopies=nmatch)
    matchbin.xbin = hiibin[m1].xbin
    for ii = 0, nmatch-1 do matchbin[ii].medy = $
      im_weighted_mean([hiibin[m1[ii]].medy,agnbin[m1[ii]].medy],$
      weights=[0.8,0.2])
;   struct_print, matchbin
    
; make some plots
    ohrange = [8.2,9.4]
    massrange = [8.2,11.9]
    residrange = [-0.8,0.4]
    drange = [-0.05,1.4]
;   levels = [0.01,0.05,0.25,0.5,0.75,0.95,0.99]
    levels = [0.5,0.75,0.9,0.975]
    
    psfile = pspath+'mzerror_agncontam'+suffix
    im_plotconfig, 1, pos, psfile=psfile, xspace=0.2, $
      xmargin=[1.1,1.2], width=[4.0,4.0], height=4.0
; MZ relation
    mzsdss_hogg_scatterplot, hiimass.mass_avg, hiioh.log12oh_t04, $
      position=pos[*,0], xsty=1, ysty=1, xrange=massrange, yrange=ohrange, $
      xtitle=mz_masstitle(), ytitle=mz_ohtitle(), levels=levels
;   djs_oplot, mzfit.bin_mass, mzfit.bin_oh, psym=6, $
;     symsize=1.3, color='blue', thick=4
    djs_oplot, mzfit.bin_mass, mzfit.bin_oh, $
      line=0, thick=10, color=fsc_color('orange red',10)
;   djs_oplot, mzfit.bin_mass, poly(mzfit.bin_mass,mzfit.coeff), $
;     line=0, thick=10, color=fsc_color('goldenrod',10)

    mzsdss_hogg_scatterplot, agnmass.mass_avg, agnoh.log12oh_t04, $
      position=pos[*,0], /noerase, xsty=5, ysty=5, $
      xrange=massrange, yrange=ohrange, /nooutlier, $
      /nogrey, ccolor=fsc_color('sienna',11), cthick=4, levels=levels
; money plot
    mzsdss_hogg_scatterplot, hiiispec[hiiindx].bpt_d, hiidist, $
      position=pos[*,1], /noerase, xsty=1, ysty=1, $
      xrange=drange, yrange=residrange, $
      xtitle='D (dex)', ytitle='', ytickname=replicate(' ',10), levels=levels
    axis, /yaxis, ysty=1, ytitle='', yrange=residrange
    xyouts, pos[2,1]+0.08, (pos[3,1]-pos[1,1])/2.0+pos[1,1], $
      'log (O/H) Residuals', /norm, align=0.5, orientation=270

    mzsdss_hogg_scatterplot, agnispec[agnindx].bpt_d, agndist, $
      position=pos[*,1], /noerase, xsty=5, ysty=5, $
      xrange=drange, yrange=residrange, /nooutlier, $
      /nogrey, ccolor=fsc_color('sienna',11), cthick=4, levels=levels

    djs_oplot, hiibin.xbin, hiibin.medy, line=0, $
      color=fsc_color('orange red',10), thick=8
    djs_oplot, agnbin.xbin, agnbin.medy, line=5, $
      color='navy', thick=10
    djs_oplot, matchbin.xbin, matchbin.medy, line=3, $
      color='red', thick=10
;   djs_oplot, hiiispec.bpt_d, hiidist, psym=3, color='red'
    im_plotconfig, /psclose

stop
return
end
    
