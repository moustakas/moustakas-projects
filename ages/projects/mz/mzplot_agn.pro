pro mzplot_agn, ps=ps
; jm09apr26nyu - AGN plots

    mzpath = mz_path()
    qapath = mzpath+'qaplots/'
    if keyword_set(ps) then begin
       pspath = qapath
       suffix = '.ps'
    endif else begin
       pspath = mz_path(/paper)
       suffix = '.eps'
    endelse

    oiihbaxis = findgen(((1.449)-(-1.5))/0.01)*0.01+(-1.5)
    oiiihblama = 0.14/(oiihbaxis-1.45)+0.83

    r23axis = findgen(((1.5)-(-0.1))/0.01)*0.01+(-0.1)
    o32lama = 1.5/(r23axis-1.7)+2.4
    ocor = 1.0+1.0/2.984

; ------------------------------------------------------------
; Yan diagnostic diagram - AGES & SDSS
    ages = read_mz_sample(/mz_ispec)
    agesanc = read_mz_sample(/mz_ancillary)
    sdss = read_mz_sample(/mz_ispec,/sdss)
    sdssanc = read_mz_sample(/mz_ancillary,/sdss)

    xrange = [0.15,1.75]
    yrange = [-1.0,1.2]
    xtitle = 'U - B'
    ytitle = textoidl('log ([O III] \lambda5007/H\beta)')
    
    sfpsym = 15  & sfcolor = 'dodger blue'  & sfsymsize = 0.5
    agnpsym = 15 & agncolor = 'dark grey'   & agnsymsize = 0.5
    unkpsym = 16 & unkcolor = 'navy' & unksymsize = 0.4
    otherpsym = 6 & othercolor = 'tan' & unksymsize = 0.6

    ages_outsymsize = 0.4
    sdss_outsymsize = 0.2
    ages_levels = [0.1,0.25,0.5,0.75,0.95]
    sdss_levels = [0.1,0.25,0.5,0.75,0.95]

    ubaxis = range(0.0,2.0,500)
    yan = (1.4-1.2*ubaxis)>(-0.1)
    
    psfile = pspath+'ages_sdss_yan'+suffix
    im_plotconfig, 1, pos, psfile=psfile, xmargin=[1.3,0.3]

; SDSS     
    ub = sdssanc.k_ubvrijhk_absmag_00[0]-sdssanc.k_ubvrijhk_absmag_00[1]
    agn = where(sdss.bpt_agn eq 1 and sdss.bpt_oiii_hb gt -900,nagn)
    sf = where(sdss.bpt_agn eq 0 and sdss.bpt_oiii_hb gt -900,nsf)
    unk = where(sdss.bpt_agn eq -1 and sdss.bpt_oiii_hb gt -900,nunk)

    mzplot_scatterplot, ub[sf], sdss[sf].bpt_oiii_hb, position=pos[*,0], $
      xsty=1, ysty=1, xtitle=xtitle, ytitle=ytitle, $
      levels=sdss_levels, cthick=6, /noout, $
      xrange=xrange, yrange=yrange, /sdss, /nogrey, ccolor=fsc_color(sfcolor,101), $
      outpsym=symcat(sfpsym,thick=8), outcolor=fsc_color(sfcolor,101), outsymsize=sdss_outsymsize
    mzplot_scatterplot, ub[agn], sdss[agn].bpt_oiii_hb, position=pos[*,0], $
      xsty=9, ysty=9, xtitle='', ytitle='', /overplot, $
      levels=sdss_levels, cthick=6, /noout, cline=0, $
      xrange=xrange, yrange=yrange, /sdss, /nogrey, ccolor=fsc_color(agncolor,101), $
      outpsym=symcat(agnpsym,thick=8), outcolor=fsc_color(agncolor,101), outsymsize=sdss_outsymsize
; for clarity do not plot the unclassified objects - 
;   djs_oplot, ub[unk], sdss[unk].bpt_oiii_hb, psym=symcat(unkpsym,thick=4), $
;     symsize=unksymsize, color=fsc_color(unkcolor,100)
    djs_oplot, ubaxis, yan, line=0, thick=7
    legend, 'SDSS', /left, /bottom, box=0, margin=0, charsize=1.6
;   legend, 'SDSS - 0.03<z<0.25', /right, /top, box=0, margin=0, charsize=1.6
    
    xyouts, 0.5, 0.9, 'AGN', orientation=-40, align=0.5, charsize=1.3
    xyouts, 0.4, 0.75, 'Star-Forming', orientation=-40, align=0.5, charsize=1.3

; AGES    
    ub = agesanc.k_ubvrijhk_absmag_00[0]-agesanc.k_ubvrijhk_absmag_00[1]
    agn = where(ages.bpt_agn eq 1 and ages.bpt_oiii_hb gt -900,nagn)
    sf = where(ages.bpt_agn eq 0 and ages.bpt_oiii_hb gt -900,nsf)
    unk = where(ages.bpt_agn eq -1 and ages.bpt_oiii_hb gt -900,nunk)
;   other = where((ages.bpt_agn eq -1) and (ages.radio_agn eq 1))
;   other = where((ages.bpt_agn eq -1) and ($
;     (ages.irac_agn eq 1) or (ages.radio_agn eq 1)))
    other = where((ages.bpt_agn eq -1) and ((ages.xray_agn eq 1) or $
      (ages.irac_agn eq 1) or (ages.radio_agn eq 1)) and (ages.bpt_oiii_hb gt -900.0))

    mzplot_scatterplot, ub[sf], ages[sf].bpt_oiii_hb, /noerase, position=pos[*,1], $
      xsty=1, ysty=1, xtitle=xtitle, ytitle='', ytickname=replicate(' ',10), $
      levels=ages_levels, cthick=6, /noout, $
      xrange=xrange, yrange=yrange, /nogrey, ccolor=fsc_color(sfcolor,101), $
      outpsym=symcat(sfpsym,thick=8), outcolor=fsc_color(sfcolor,101), outsymsize=ages_outsymsize
    mzplot_scatterplot, ub[agn], ages[agn].bpt_oiii_hb, position=pos[*,0], $
      xsty=9, ysty=9, xtitle='', ytitle='', /overplot, $
      levels=ages_levels, cthick=6, /noout, cline=0, $
      xrange=xrange, yrange=yrange, /nogrey, ccolor=fsc_color(agncolor,101), $
      outpsym=symcat(agnpsym,thick=8), outcolor=fsc_color(agncolor,101), outsymsize=ages_outsymsize
    djs_oplot, ub[unk], ages[unk].bpt_oiii_hb, psym=symcat(unkpsym,thick=4), $
      symsize=unksymsize, color=fsc_color(unkcolor,100)
    djs_oplot, ub[other], ages[other].bpt_oiii_hb, psym=symcat(otherpsym,thick=4), $
      symsize=othersymsize, color=fsc_color(othercolor,100)
    
;   djs_plot, [0], [0], /nodata, /noerase, position=pos[*,1], xsty=1, ysty=1, $
;     xtitle=xtitle, ytitle='', ytickname=replicate(' ',10), $
;     xrange=xrange, yrange=yrange
;   djs_oplot, ub[sf], ages[sf].bpt_oiii_hb, psym=symcat(sfpsym,thick=4), $
;     symsize=sfsymsize, color=fsc_color(sfcolor,100)
;   djs_oplot, ub[agn], ages[agn].bpt_oiii_hb, psym=symcat(agnpsym,thick=4), $
;     symsize=agnsymsize, color=fsc_color(agncolor,100)
;   djs_oplot, ub[unk], ages[unk].bpt_oiii_hb, psym=symcat(unkpsym,thick=4), $
;     symsize=unksymsize, color=fsc_color(unkcolor,100)

    djs_oplot, ubaxis, yan, line=0, thick=7
    legend, 'AGES', /left, /bottom, box=0, margin=0, charsize=1.6
;   legend, 'AGES - 0.05<z<0.75', /right, /top, box=0, margin=0, charsize=1.6
;   xyouts, 1.55, 1.1, 'AGN', align=0.5, charsize=1.6
;   xyouts, 0.4, -0.6, 'Star-!cForming', align=0.5, charsize=1.6
    
    im_plotconfig, /psclose, psfile=psfile, gzip=keyword_set(ps)

stop    
    
; ------------------------------------------------------------
; [N II]/Ha versus [O III]/Hb - AGES
    ages = read_mz_sample(/mz_ispec)
    agesanc = read_mz_sample(/mz_ancillary)

    xrange = [-1.8,0.6]
    yrange = [-1.4,1.1]
    xtitle = textoidl('log ([N II] \lambda6584/H\alpha)')
    ytitle = textoidl('log ([O III] \lambda5007/H\beta)')

    psfile = pspath+'ages_bpt'+suffix
    im_plotconfig, 0, pos, psfile=psfile, charsize=2.0, xmargin=[1.6,0.4], $
      height=6

    good = where(ages.bpt_agn ne -1)
    mzplot_scatterplot, ages[good].bpt_nii_ha, ages[good].bpt_oiii_hb, $
      position=pos, xsty=1, ysty=1, xtitle=xtitle, ytitle=ytitle, $
      xrange=xrange, yrange=yrange, npix=24

; case 1    
    case1 = where((ages.bpt_nii_ha_limit gt -900.0) and $ ; [NII] limit, [OIII] limit
      (ages.bpt_oiii_hb_limit gt -900.0),ncase1)                                             

    plotsym, 6, 1.5, color=djs_icolor('orange'), thick=3
    djs_oplot, ages[case1].bpt_nii_ha_limit, ages[case1].bpt_oiii_hb_limit, psym=8
    plotsym, 1, 1.5, color=djs_icolor('orange'), thick=3
    djs_oplot, ages[case1].bpt_nii_ha_limit, ages[case1].bpt_oiii_hb_limit, psym=8
    
; case 2
    case2 = where((ages.bpt_nii_ha_limit gt -900.0) and $ ; [NII] limit, [OIII] good
      (ages.bpt_oiii_hb gt -900.0),ncase2)                                               

    plotsym, 6, 1.5, color=djs_icolor('blue'), thick=3
    djs_oplot, ages[case2].bpt_nii_ha_limit, ages[case2].bpt_oiii_hb, psym=8

; case 3    
    case3 = where((ages.bpt_nii_ha gt -900.0) and $ ; [NII] good, [OIII] limit
      (ages.bpt_oiii_hb_limit gt -900.0),ncase3) 

    plotsym, 1, 1.5, color=djs_icolor('dark green'), thick=3
    djs_oplot, ages[case3].bpt_nii_ha, ages[case3].bpt_oiii_hb_limit, psym=8
    
    models = kewley_bpt_lines(/kauffmann,_extra=extra)
    oplot, models.x_nii, models.y_nii, line=0, thick=6.0
    djs_oplot, -0.3*[1,1], !y.crange, line=1, thick=6
    
    xyouts, 0.4, 0.45, 'AGN', align=0.5, charsize=1.8
    xyouts, -1.4, -0.55, 'Star-Forming', align=0.5, charsize=1.8
    
    im_plotconfig, /psclose, psfile=psfile, gzip=keyword_set(ps)

; ------------------------------------------------------------
; EW([O III])/EW(Hb) vs EW([O II])/EW(Hb) - AGES + SDSS
    
    psfile = qapath+'ages_sdss_oiiihb_vs_oiihb.ps'
    im_plotconfig, 1, pos, psfile=psfile, charsize=1.5

; AGES    
    unkages = where((strmatch(ages.class,'*UNKNOWN*') eq 1B),nunkages)
    agnages = where((strmatch(ages.class,'*UNKNOWN*') eq 0B) and $
      (strmatch(ages.class,'*BPT_AGN*') eq 1B),nagnages)
    sfages = where((strmatch(ages.class,'*UNKNOWN*') eq 0B) and $
      (strmatch(ages.class,'*BPT_AGN*') eq 0B),nsfages)
    help, nsfages, nagnages, nunkages

    xsfages = alog10(ages[sfages].oii_3727_ew[0]/ages[sfages].h_beta_ew[0])
    ysfages = alog10(ages[sfages].oiii_5007_ew[0]/ages[sfages].h_beta_ew[0])
    xagnages = alog10(ages[agnages].oii_3727_ew[0]/ages[agnages].h_beta_ew[0])
    yagnages = alog10(ages[agnages].oiii_5007_ew[0]/ages[agnages].h_beta_ew[0])
    xunkages = alog10(ages[unkages].oii_3727_ew[0]/ages[unkages].h_beta_ew[0])
    yunkages = alog10(ages[unkages].oiii_5007_ew[0]/ages[unkages].h_beta_ew[0])

; SDSS    
    unksdss = where((strmatch(sdss.class,'*UNKNOWN*') eq 1B),nunksdss)
    agnsdss = where((strmatch(sdss.class,'*UNKNOWN*') eq 0B) and $
      (strmatch(sdss.class,'*BPT_AGN*') eq 1B),nagnsdss)
    sfsdss = where((strmatch(sdss.class,'*UNKNOWN*') eq 0B) and $
      (strmatch(sdss.class,'*BPT_AGN*') eq 0B),nsfsdss)
    help, nsfsdss, nagnsdss, nunksdss

    xsfsdss = alog10(sdss[sfsdss].oii_3727_ew[0]/sdss[sfsdss].h_beta_ew[0])
    ysfsdss = alog10(sdss[sfsdss].oiii_5007_ew[0]/sdss[sfsdss].h_beta_ew[0])
    xagnsdss = alog10(sdss[agnsdss].oii_3727_ew[0]/sdss[agnsdss].h_beta_ew[0])
    yagnsdss = alog10(sdss[agnsdss].oiii_5007_ew[0]/sdss[agnsdss].h_beta_ew[0])
;   xunksdss = alog10(sdss[unksdss].oii_3727_ew[0]/sdss[unksdss].h_beta_ew[0])
;   yunksdss = alog10(sdss[unksdss].oiii_5007_ew[0]/sdss[unksdss].h_beta_ew[0])

; what fraction of the AGN/SF are classified as AGN using this method? 
    splog, '## AGES ##'
    oiiihblama_agn = 0.14/(xagnages-1.45)+0.83
    lama_agn = where(yagnages gt oiiihblama_agn,nlama_agn)
    oiiihblama_sf = 0.14/(xsfages-1.45)+0.83
    lama_sf = where(ysfages gt oiiihblama_sf,nlama_sf)
    oiiihblama_unk = 0.14/(xunkages-1.45)+0.83
    lama_unk = where(yunkages gt oiiihblama_unk,nlama_unk)

    splog, 'Fraction of AGN above curve = ', nlama_agn, nagnages, 100.0*nlama_agn/float(nagnages)
    splog, 'Fraction of SF above curve = ', nlama_sf, nsfages, 100.0*nlama_sf/float(nsfages)
    splog, 'Fraction of Unclassified above curve = ', nlama_unk, nunkages, 100.0*nlama_unk/float(nunkages)
    
    splog, '## SDSS ##'
    oiiihblama_agn = 0.14/(xagnsdss-1.45)+0.83
    lama_agn = where(yagnsdss gt oiiihblama_agn,nlama_agn)
    oiiihblama_sf = 0.14/(xsfsdss-1.45)+0.83
    lama_sf = where(ysfsdss gt oiiihblama_sf,nlama_sf)
;   oiiihblama_unk = 0.14/(xunksdss-1.45)+0.83
;   lama_unk = where(yunksdss gt oiiihblama_unk,nlama_unk)

    splog, 'Fraction of AGN above curve = ', nlama_agn, nagnsdss, 100.0*nlama_agn/float(nagnsdss)
    splog, 'Fraction of SF above curve = ', nlama_sf, nsfsdss, 100.0*nlama_sf/float(nsfsdss)
;   splog, 'Fraction of Unclassified above curve = ', nlama_unk, nunksdss, 100.0*nlama_unk/float(nunksdss)
    
; now make the plot    
    xtitle = textoidl('log {EW([O II] \lambda3727)/EW(H\beta)}')
    ytitle = textoidl('log {EW([O III] \lambda5007)/EW(H\beta)}')

    xrange = [-0.4,1.6] ; oiihbrange
    yrange = [-1.1,1.4] ; oiiihbrange

; SDSS    
    mzsdss_hogg_scatterplot, xsfsdss, ysfsdss, position=pos[*,0], $
      xsty=1, ysty=1, xtitle=xtitle, ytitle=ytitle, $
      xrange=xrange, yrange=yrange, /nooutliers
    these = long(randomu(seed,fix(0.2*nagnsdss))*nagnsdss) ; plot 20%
    agnoplot, xagnsdss[these], yagnsdss[these], /sdss
; don't plot the unclassified SDSS galaxies - we toss them out
;   unkoplot, xunksdss, yunksdss, /sdss

    djs_oplot, oiihbaxis, oiiihblama, line=0, thick=6.0
;   djs_oplot, oiihbaxis+0.1, oiiihblama+0.1, line=5, thick=6.0
;   djs_oplot, oiihbaxis-0.1, oiiihblama-0.1, line=5, thick=6.0

    im_legend, 'SDSS - 0.03<z<0.25', /left, /top, box=0, margin=0
    
; label the regions    
    xyouts, 1.32, 0.95, 'AGN', align=0.5
    xyouts, 0.9, -0.9, 'Star-Forming', align=0.5
    
; AGES
    mzages_hogg_scatterplot, xsfages, ysfages, /noerase, position=pos[*,1], $
      xsty=1, ysty=1, xtitle=xtitle, ytitle='', ytickname=replicate(' ',10), $
      xrange=xrange, yrange=yrange, /nooutliers
;   djs_oplot, !x.crange, [0.75,0.75], line=5, thick=6.0
    agnoplot, xagnages, yagnages
;   unkoplot, xunkages, yunkages
    these = long(randomu(seed,fix(0.2*nunkages))*nunkages) ; plot 20%
    unkoplot, xunkages[these], yunkages[these]

    djs_oplot, oiihbaxis, oiiihblama, line=0, thick=6.0
;   djs_oplot, oiihbaxis+0.1, oiiihblama+0.1, line=5, thick=6.0
;   djs_oplot, oiihbaxis-0.1, oiiihblama-0.1, line=5, thick=6.0

    im_legend, 'AGES - 0.05<z<0.75', /left, /top, box=0, margin=0
    
; label the regions    
    xyouts, 1.32, 0.95, 'AGN', align=0.5
    xyouts, 0.9, -0.9, 'Star-Forming', align=0.5
    
    im_plotconfig, /psclose, psfile=psfile, /gzip

; ------------------------------------------------------------
; EW(R23) vs EW(O32) - AGES + SDSS
    
    psfile = qapath+'ages_sdss_r23_vs_o32.ps'
    im_plotconfig, 1, pos, psfile=psfile, charsize=1.5

; AGES    
    unkages = where((strmatch(ages.class,'*UNKNOWN*') eq 1B),nunkages)
    agnages = where((strmatch(ages.class,'*UNKNOWN*') eq 0B) and $
      (strmatch(ages.class,'*BPT_AGN*') eq 1B),nagnages)
    sfages = where((strmatch(ages.class,'*UNKNOWN*') eq 0B) and $
      (strmatch(ages.class,'*BPT_AGN*') eq 0B),nsfages)
    help, nsfages, nagnages, nunkages

    xsfages = alog10((ages[sfages].oii_3727_ew[0]+ocor*ages[sfages].oiii_5007_ew[0])/ages[sfages].h_beta_ew[0])
    ysfages = alog10((ocor*ages[sfages].oiii_5007_ew[0])/ages[sfages].oii_3727_ew[0])
    xagnages = alog10((ages[agnages].oii_3727_ew[0]+ocor*ages[agnages].oiii_5007_ew[0])/ages[agnages].h_beta_ew[0])
    yagnages = alog10((ocor*ages[agnages].oiii_5007_ew[0])/ages[agnages].oii_3727_ew[0])
    xunkages = alog10((ages[unkages].oii_3727_ew[0]+ocor*ages[unkages].oiii_5007_ew[0])/ages[unkages].h_beta_ew[0])
    yunkages = alog10((ocor*ages[unkages].oiii_5007_ew[0])/ages[unkages].oii_3727_ew[0])

; SDSS    
    unksdss = where((strmatch(sdss.class,'*UNKNOWN*') eq 1B),nunksdss)
    agnsdss = where((strmatch(sdss.class,'*UNKNOWN*') eq 0B) and $
      (strmatch(sdss.class,'*BPT_AGN*') eq 1B),nagnsdss)
    sfsdss = where((strmatch(sdss.class,'*UNKNOWN*') eq 0B) and $
      (strmatch(sdss.class,'*BPT_AGN*') eq 0B),nsfsdss)
    help, nsfsdss, nagnsdss, nunksdss

    xsfsdss = alog10((sdss[sfsdss].oii_3727_ew[0]+ocor*sdss[sfsdss].oiii_5007_ew[0])/sdss[sfsdss].h_beta_ew[0])
    ysfsdss = alog10((ocor*sdss[sfsdss].oiii_5007_ew[0])/sdss[sfsdss].oii_3727_ew[0])
    xagnsdss = alog10((sdss[agnsdss].oii_3727_ew[0]+ocor*sdss[agnsdss].oiii_5007_ew[0])/sdss[agnsdss].h_beta_ew[0])
    yagnsdss = alog10((ocor*sdss[agnsdss].oiii_5007_ew[0])/sdss[agnsdss].oii_3727_ew[0])
;   xunksdss = alog10((sdss[unksdss].oii_3727_ew[0]+ocor*sdss[unksdss].oiii_5007_ew[0])/sdss[unksdss].h_beta_ew[0])
;   yunksdss = alog10((ocor*sdss[unksdss].oiii_5007_ew[0])/sdss[unksdss].oii_3727_ew[0])

; what fraction of the AGN/SF are classified as AGN using this method? 
    splog, '## AGES ##'
    o32lama_agn = 1.5/(xagnages-1.7)+2.4
    lama_agn = where(yagnages gt o32lama_agn,nlama_agn)
    o32lama_sf = 1.5/(xsfages-1.7)+2.4
    lama_sf = where(ysfages gt o32lama_sf,nlama_sf)
    o32lama_unk = 1.5/(xunkages-1.7)+2.4
    lama_unk = where(yunkages gt o32lama_unk,nlama_unk)

    splog, 'Fraction of AGN above curve = ', nlama_agn, nagnages, 100.0*nlama_agn/float(nagnages)
    splog, 'Fraction of SF above curve = ', nlama_sf, nsfages, 100.0*nlama_sf/float(nsfages)
    splog, 'Fraction of Unclassified above curve = ', nlama_unk, nunkages, 100.0*nlama_unk/float(nunkages)
    
    splog, '## SDSS ##'
    o32lama_agn = 1.5/(xagnsdss-1.7)+2.4
    lama_agn = where(yagnsdss gt o32lama_agn,nlama_agn)
    o32lama_sf = 1.5/(xsfsdss-1.7)+2.4
    lama_sf = where(ysfsdss gt o32lama_sf,nlama_sf)
;   o32lama_unk = 1.5/(xunksdss-1.7)+2.4
;   lama_unk = where(yunksdss gt o32lama_unk,nlama_unk)

    splog, 'Fraction of AGN above curve = ', nlama_agn, nagnsdss, 100.0*nlama_agn/float(nagnsdss)
    splog, 'Fraction of SF above curve = ', nlama_sf, nsfsdss, 100.0*nlama_sf/float(nsfsdss)
;   splog, 'Fraction of Unclassified above curve = ', nlama_unk, nunksdss, 100.0*nlama_unk/float(nunksdss)
    
; now make the plot
    xtitle = textoidl('log EW(R_{23})')
    ytitle = textoidl('log EW(O_{32})')

    xrange = [-0.1,1.6]
    yrange = [-1.5,1.0]

; AGES
    mzages_hogg_scatterplot, xsfages, ysfages, position=pos[*,0], $
      xtitle=xtitle, ytitle=ytitle, xsty=1, ysty=1, $
      xrange=xrange, yrange=yrange, /nooutliers
    agnoplot, xagnages, yagnages
    unkoplot, xunkages, yunkages
    djs_oplot, r23axis, o32lama, line=0, thick=6.0
;   djs_oplot, r23axis+0.1, o32lama+0.1, line=5, thick=6.0
;   djs_oplot, r23axis-0.1, o32lama-0.1, line=5, thick=6.0

    im_legend, 'AGES - 0.05<z<0.75', /left, /top, box=0, charsize=1.5, margin=0

; SDSS    
    mzsdss_hogg_scatterplot, xsfsdss, ysfsdss, /noerase, position=pos[*,1], $
      xsty=1, ysty=1, xtitle=xtitle, ytitle='', ytickname=replicate(' ',10), $
      xrange=xrange, yrange=yrange, /nooutliers
    these = long(randomu(seed,fix(0.2*nagnsdss))*nagnsdss) ; plot 20%
    agnoplot, xagnsdss[these], yagnsdss[these], /sdss
; don't plot the unclassified SDSS galaxies - we toss them out
;   unkoplot, xunksdss, yunksdss, /sdss

    djs_oplot, r23axis, o32lama, line=0, thick=6.0
;   djs_oplot, r23axis+0.1, o32lama+0.1, line=5, thick=6.0
;   djs_oplot, r23axis-0.1, o32lama-0.1, line=5, thick=6.0

    im_legend, 'SDSS - 0.03<z<0.25', /left, /top, box=0, charsize=1.5, margin=0

    im_plotconfig, /psclose, psfile=psfile, /gzip

return
end
    
;
;
;; ------------------------------------------------------------
;; [N II]/Ha versus [O III]/Hb - AGES
;    psfile = pspath+'ages_bpt'+suffix
;    im_plotconfig, 0, pos, psfile=psfile, charsize=2.0, xmargin=[1.6,0.4]
;
;    unkages = where(strmatch(ages.class,'*UNKNOWN*'),$
;      nunkages,comp=good,ncomp=ngood)
;    ages1 = ages[good]
;    nages = n_elements(ages)
;
;    agnages = where(strmatch(ages1.class,'*BPT_AGN*'),nagnages)
;    sfages = where((strmatch(ages1.class,'*BPT_AGN*') eq 0B),nsfages)
;    help, nages, nunkages, ngood, nsfages, nagnages
;
;    agnlama = where(strmatch(ages1.class,'*LAMA_AGN*'))
;    agnxray = where(strmatch(ages1.class,'*XRAY_AGN*'))
;    agnnev = where(strmatch(ages1.class,'*NEV_AGN*'))
;    agnirac = where(strmatch(ages1.class,'*IRAC_AGN*'))
;    agnassef = where(strmatch(ages1.class,'*ASSEF_AGN*'))
;
;    xages = alog10(ages1.nii_6584[0]/ages1.h_alpha[0])
;    yages = alog10(ages1.oiii_5007[0]/ages1.h_beta[0])
;
;    xsfages = alog10(ages1[sfages].nii_6584[0]/ages1[sfages].h_alpha[0])
;    ysfages = alog10(ages1[sfages].oiii_5007[0]/ages1[sfages].h_beta[0])
;
;    xagnages = alog10(ages1[agnages].nii_6584[0]/ages1[agnages].h_alpha[0])
;    yagnages = alog10(ages1[agnages].oiii_5007[0]/ages1[agnages].h_beta[0])
;    
;    xtitle = textoidl('log ([N II] \lambda6584/H\alpha)')
;    ytitle = textoidl('log ([O III] \lambda5007/H\beta)')
;
;    xrange = [-1.8,0.7]
;    yrange = [-1.1,1.3]
;
;    djs_plot, [0], [0], /nodata, position=pos[*,0], xsty=1, ysty=1, $
;      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange
;    sfoplot, xsfages, ysfages
;    agnoplot, xagnages, yagnages
;
;    djs_oplot, xages[agnxray], yages[agnxray], psym=symcat(6,thick=3), $
;      symsize=0.9, color='blue'
;    djs_oplot, xages[agnirac], yages[agnirac], psym=symcat(6,thick=3), $
;      symsize=0.9, color='yellow'
;    djs_oplot, xages[agnlama], yages[agnlama], psym=symcat(6,thick=3), $
;      symsize=0.9, color='orange'
;;   djs_oplot, xages[agnassef], yages[agnassef], psym=symcat(6,thick=3), $
;;     symsize=0.9, color='red'
;    djs_oplot, xages[agnnev], yages[agnnev], psym=symcat(6,thick=3), $
;      symsize=0.9, color='green'
;    
;; overplot the Kauffmann and Kewley mixing lines
;    models = kewley_bpt_lines(/kauffmann,_extra=extra)
;    oplot, models.x_nii, models.y_nii, line=0, thick=6.0
;
;    models = kewley_bpt_lines(_extra=extra)
;    oplot, models.x_nii, models.y_nii, line=5, thick=6.0
;
;; label the regions
;    im_legend, 'AGES - 0.05<z<0.4', /left, /top, box=0, charsize=1.8
;    xyouts, 0.45, 0.45, 'AGN', align=0.5, charsize=1.8
;    xyouts, -1.25, -0.4, 'Star-Forming', align=0.5, charsize=1.8
;    
;    im_plotconfig, /psclose
;
;stop    
;    


;; ------------------------------------------------------------
;; Yan diagnostic diagram - AGES & SDSS
;    xrange = [-1.8,0.7]
;    yrange = [-1.1,1.3]
;
;    sfpsym = 9   & sfcolor = 'dodger blue'  & sfsymsize = 0.5
;    agnpsym = 15 & agncolor = 'firebrick'   & agnsymsize = 0.5
;    unkpsym = 14 & agncolor = 'royal blue ' & agnsymsize = 0.5
;    
;    psfile = pspath+'ages_bpt'+suffix
;    im_plotconfig, 0, pos, psfile=psfile, charsize=2.0, xmargin=[1.6,0.4]
;
;    agn = where(ages.bpt_agn eq 1)
;    sf = where(ages.bpt_agn eq 0)
;    unk = where(ages.bpt_agn eq -1)
;    
;    djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
;      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange
;    djs_oplot, ages[agn].bpt_nii_ha, agn[agn].bpt_oiii_hb, $
;      psym=symcat(agnpsym,thick=4), symsize=agnsymsize, $
;      color=fsc_color(agncolor,100)
;    djs_oplot, ages[sf].bpt_nii_ha, sf[sf].bpt_oiii_hb, $
;      psym=symcat(sfpsym,thick=4), symsize=sfsymsize, $
;      color=fsc_color(sfcolor,100)
;    
;    
;    unkages = where(strmatch(ages.class,'*UNKNOWN*'),$
;      nunkages,comp=good,ncomp=ngood)
;    ages1 = ages[good]
;    nages = n_elements(ages)
;
;    agnages = where(strmatch(ages1.class,'*BPT_AGN*'),nagnages)
;    sfages = where((strmatch(ages1.class,'*BPT_AGN*') eq 0B),nsfages)
;    help, nages, nunkages, ngood, nsfages, nagnages
;
;    agnlama = where(strmatch(ages1.class,'*LAMA_AGN*'))
;    agnxray = where(strmatch(ages1.class,'*XRAY_AGN*'))
;    agnnev = where(strmatch(ages1.class,'*NEV_AGN*'))
;    agnirac = where(strmatch(ages1.class,'*IRAC_AGN*'))
;    agnassef = where(strmatch(ages1.class,'*ASSEF_AGN*'))
;
;    xages = alog10(ages1.nii_6584[0]/ages1.h_alpha[0])
;    yages = alog10(ages1.oiii_5007[0]/ages1.h_beta[0])
;
;    xsfages = alog10(ages1[sfages].nii_6584[0]/ages1[sfages].h_alpha[0])
;    ysfages = alog10(ages1[sfages].oiii_5007[0]/ages1[sfages].h_beta[0])
;
;    xagnages = alog10(ages1[agnages].nii_6584[0]/ages1[agnages].h_alpha[0])
;    yagnages = alog10(ages1[agnages].oiii_5007[0]/ages1[agnages].h_beta[0])
;    
;    xtitle = textoidl('log ([N II] \lambda6584/H\alpha)')
;    ytitle = textoidl('log ([O III] \lambda5007/H\beta)')
;
;    sfoplot, xsfages, ysfages
;    agnoplot, xagnages, yagnages
;
;    djs_oplot, xages[agnxray], yages[agnxray], psym=symcat(6,thick=3), $
;      symsize=0.9, color='blue'
;    djs_oplot, xages[agnirac], yages[agnirac], psym=symcat(6,thick=3), $
;      symsize=0.9, color='yellow'
;    djs_oplot, xages[agnlama], yages[agnlama], psym=symcat(6,thick=3), $
;      symsize=0.9, color='orange'
;;   djs_oplot, xages[agnassef], yages[agnassef], psym=symcat(6,thick=3), $
;;     symsize=0.9, color='red'
;    djs_oplot, xages[agnnev], yages[agnnev], psym=symcat(6,thick=3), $
;      symsize=0.9, color='green'
;    
;; overplot the Kauffmann and Kewley mixing lines
;    models = kewley_bpt_lines(/kauffmann,_extra=extra)
;    oplot, models.x_nii, models.y_nii, line=0, thick=6.0
;
;    models = kewley_bpt_lines(_extra=extra)
;    oplot, models.x_nii, models.y_nii, line=5, thick=6.0
;
;; label the regions
;    im_legend, 'AGES - 0.05<z<0.4', /left, /top, box=0, charsize=1.8
;    xyouts, 0.45, 0.45, 'AGN', align=0.5, charsize=1.8
;    xyouts, -1.25, -0.4, 'Star-Forming', align=0.5, charsize=1.8
;    
;    im_plotconfig, /psclose
;
;stop    
;    
;;; ------------------------------------------------------------
;;; [N II]/Ha versus [O III]/Hb - SDSS
;;    
;;    psfile = pspath+'sdss_bpt'+suffix
;;    im_plotconfig, 0, pos, psfile=psfile, xmargin=[1.3,0.4], $
;;      width=6.8, height=6.8
;;
;;    xsdss = alog10(sdss.nii_6584[0]/sdss.h_alpha[0])
;;    ysdss = alog10(sdss.oiii_5007[0]/sdss.h_beta[0])
;;
;;;    unksdss = where(strmatch(sdss.class,'*UNKNOWN*'),$
;;;      nunksdss,comp=good,ncomp=ngood)
;;;    sdss1 = sdss[good]
;;;    nsdss = n_elements(sdss)
;;;    agnsdss = where(strmatch(sdss1.class,'*BPT_AGN*'),nagnsdss)
;;;    sfsdss = where((strmatch(sdss1.class,'*BPT_AGN*') eq 0B),nsfsdss)
;;;    help, nsdss, nunksdss, ngood, nsfsdss, nagnsdss
;;;
;;;    xsdss = alog10(sdss1.nii_6584[0]/sdss1.h_alpha[0])
;;;    ysdss = alog10(sdss1.oiii_5007[0]/sdss1.h_beta[0])
;;;
;;;    xsfsdss = alog10(sdss1[sfsdss].nii_6584[0]/sdss1[sfsdss].h_alpha[0])
;;;    ysfsdss = alog10(sdss1[sfsdss].oiii_5007[0]/sdss1[sfsdss].h_beta[0])
;;;    xagnsdss = alog10(sdss1[agnsdss].nii_6584[0]/sdss1[agnsdss].h_alpha[0])
;;;    yagnsdss = alog10(sdss1[agnsdss].oiii_5007[0]/sdss1[agnsdss].h_beta[0])
;;    
;;    xtitle = textoidl('log ([N II] \lambda6584/H\alpha)')
;;    ytitle = textoidl('log ([O III] \lambda5007/H\beta)')
;;
;;    xrange = [-1.8,0.7]
;;    yrange = [-1.1,1.3]
;;
;;    mzsdss_hogg_scatterplot, xsdss, ysdss, position=pos[*,0], $
;;      xsty=1, ysty=1, xtitle=xtitle, ytitle=ytitle, $
;;      xrange=xrange, yrange=yrange;, /nooutliers
;;;   mzsdss_hogg_scatterplot, xagnsdss, yagnsdss, /noerase, position=pos[*,0], $
;;;     xsty=1, ysty=1, xtitle='', ytitle='', xrange=xrange, yrange=yrange, /overplot
;;;   agnoplot, xagnsdss, yagnsdss, /sdss
;;
;;; overplot the Kauffmann and Kewley mixing lines
;;    models = kewley_bpt_lines(/kauffmann,_extra=extra)
;;    oplot, models.x_nii, models.y_nii, line=0, thick=6.0
;;
;;;   models = kewley_bpt_lines(_extra=extra)
;;;   oplot, models.x_nii, models.y_nii, line=5, thick=6.0
;;
;;; label the regions
;;    im_legend, 'SDSS - 0.03<z<0.25', /left, /top, box=0, charsize=1.8
;;    xyouts, 0.45, 0.45, 'AGN', align=0.5, charsize=1.8
;;    xyouts, -1.25, -0.4, 'Star-Forming', align=0.5, charsize=1.8
;;    
;;    im_plotconfig, /psclose
;;
;;stop    
;;    
;;pro sfoplot, str, sdss=sdss, _extra=extra
;;    if keyword_set(sdss) then symsize = 0.2 else symsize = 0.5
;;    djs_oplot, str.bpt_nii_ha, str.bpt_oiii_hb, psym=symcat(9,thick=4), $
;;      symsize=symsize, color=fsc_color('dodger blue',100), _extra=extra
;;return
;;end
;;pro agnoplot, str, sdss=sdss, _extra=extra
;;    if keyword_set(sdss) then symsize = 0.2 else symsize = 0.55
;;    djs_oplot, str.bpt_nii_ha, str.bpt_oiii_hb, psym=symcat(15,thick=4), $
;;      color=fsc_color('firebrick',101), _extra=extra
;;return
;;end
;;pro unkoplot, x, y, sdss=sdss, _extra=extra
;;    if keyword_set(sdss) then symsize = 0.2 else symsize = 0.6
;;    djs_oplot, str.bpt_nii_ha, str.bpt_oiii_hb, psym=symcat(14,thick=4), $
;;      color=fsc_color('royal blue',102), _extra=extra
;;return
;;end
;;
