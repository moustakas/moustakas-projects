pro plotsings_log12oh_pmethod, ps=ps
; jm10mar09ucsd - use the SDSS sample to show that PT05 should not be
; applied to large galaxy samples
    
    if keyword_set(ps) then suffix = '.ps' else suffix = '.eps'

    metpath = sings_path(/projects)+'log12oh/'
    pspath = sings_path(/papers)+'log12oh/FIG_LOG12OH/'

; deal with the HII regions    
    allhii = read_hiiregions()
    keep = where((allhii.zstrong_r23 gt -900.0) and $
      (allhii.ewhb gt -900.0) and $
      (allhii.lit_log12oh_te gt 8.2) and $
      ((allhii.lit_t4363/allhii.lit_t4363_err gt 10.0) or $
      (allhii.lit_t5755/allhii.lit_t5755_err gt 10.0) or $
      (allhii.lit_t6312/allhii.lit_t6312_err gt 10.0)))
    hii = allhii[keep]

    texref = strtrim(hii.texref,2)
    keep = where((texref eq 'bresolin04a') or $
      (texref eq 'bresolin05a') or $
;     (texref eq 'castellanos02a') or $
      (texref eq 'garnett97a') or $
      (texref eq 'izotov04a') or $
      (texref eq 'izotov06a') or $
      (texref eq 'kenn03a') or $
      (texref eq 'kniazev04a') or $
      (texref eq 'rosolowsky08a') or $
      (texref eq 'vanzee98a'),nkeep)
    hii = hii[keep]
    splog, 'HII regions: ', nkeep

    oh = hii.lit_log12oh_te
    minoh = min(oh,mindx)
    maxoh = max(oh,maxdx)
    splog, 'Minimum O/H ', minoh, ', '+strtrim(hii[mindx].hii_galaxy,2)+', '+$
      strtrim(hii[mindx].hii_region,2)+', '+strtrim(hii[mindx].texref,2)
    splog, 'Maximum O/H ', maxoh, ', '+strtrim(hii[maxdx].hii_galaxy,2)+', '+$
      strtrim(hii[maxdx].hii_region,2)+', '+strtrim(hii[maxdx].texref,2)

; code to pick out the references I want, or to show the number of HII
; regions per reference
    texref = strtrim(hii.texref,2)
    utexref = texref[uniq(texref,sort(texref))]
    number = intarr(n_elements(utexref))
    for ii = 0, n_elements(utexref)-1 do number[ii] = $
      total(utexref[ii] eq texref)
    niceprint, number, utexref

; now the SDSS
    sdss = read_sings_log12oh_samples(/sdss)
    class = iclassification(sdss,snrcut=3.0,ratios=ratios)

    sf = where(strmatch(ratios.final_class,'*AGN*',/fold) eq 0)
;   sf = where(strtrim(ratios.final_class,2) ne 'AGN')
    sdss_nodust = iunred_linedust(sdss[sf],snrcut=3.0)
    zz = im_abundance(sdss_nodust,snrcut=3.0,$
      /justflux,/nodensity,nmonte=0)
    
; --------------------------------------------------
; R23 vs P for SDSS galaxies    
    
    xtitle = textoidl('log_{10} (R_{23})')
    ytitle = 'P'
    xrange = [-0.3,1.2]
    yrange = [0.0,1.0]
    levels = [0.5,0.75,0.95]
    cannotation = ['50%','75%','95%']
    
    psfile = pspath+'sdss_r23_vs_p'+suffix
    im_plotconfig, 0, pos, psfile=psfile
    keep = where((zz.zstrong_r23 gt -900.0) and (zz.zstrong_p gt -900.0))
    hogg_scatterplot, alog10(zz[keep].zstrong_r23), zz[keep].zstrong_p, $
      position=pos, xsty=1, ysty=1, /outliers, /internal, label=1, $
      outpsym=6, outsymsize=0.05, outcolor=djs_icolor('grey'), $
      levels=levels, cannotation=cannotation, xtitle=xtitle, ytitle=ytitle, $
      xrange=xrange, yrange=yrange
; overplot the HII regions, making the symbol size proportional to
; metallicity 
    oh = hii.lit_log12oh_te
    symsize = 1.5*sqrt(10^((oh-min(oh))/(max(oh)-min(oh))))
    djs_oplot, alog10(hii.zstrong_r23), hii.zstrong_p, psym=6, $
      color=fsc_color('firebrick',101), symsize=symsize, $
      thick=5.0
    im_legend, 'Metal-Rich HII Regions', /left, /top, box=0, $
      psym=symcat(6,thick=5.0), color='firebrick', $
      charsize=1.7

; excitation arrow    
    arrow, -0.2, 0.55, -0.2, 0.8, /data, thick=6.0, $
      hsize=-0.15, hthick=6.0
    xyouts, -0.1, 0.65, 'increasing excitation', /data, $
      orientation=90, align=0.5, charsize=1.4
; metallicity arrow    
    arrow, 1.05, 0.1, 0.7, 0.1, /data, thick=6.0, $
      hsize=-0.15, hthick=6.0
    xyouts, 0.875, 0.05, 'increasing metallicity', /data, $
      orientation=0, align=0.5, charsize=1.4
    im_plotconfig, /psclose

; --------------------------------------------------
; R23 vs EW(Hb) for SDSS galaxies 
    
    xtitle = textoidl('log_{10} (R_{23})')
    ytitle = textoidl('EW(H\beta) (\AA)')
    xrange = [-0.3,1.2]
    yrange = [0,3]
    levels = [0.5,0.75,0.95]
    cannotation = ['50%','75%','95%']
    
    psfile = pspath+'sdss_r23_vs_ewhb'+suffix
    im_plotconfig, 0, pos, psfile=psfile
    keep = where((zz.zstrong_r23 gt -900.0) and (sdss[sf].h_beta_ew[0] gt 0.0))
    hogg_scatterplot, alog10(zz[keep].zstrong_r23), alog10(sdss[sf[keep]].h_beta_ew[0]), $
      position=pos, xsty=1, ysty=1, /outliers, /internal, label=1, $
      outpsym=6, outsymsize=0.05, outcolor=djs_icolor('grey'), $
      levels=levels, cannotation=cannotation, xtitle=xtitle, ytitle=ytitle, $
      xrange=xrange, yrange=yrange
; overplot the HII regions, making the symbol size proportional to
; metallicity 
    oh = hii.lit_log12oh_te
    symsize = 1.5*sqrt(10^((oh-min(oh))/(max(oh)-min(oh))))
    djs_oplot, alog10(hii.zstrong_r23), alog10(hii.ewhb), psym=6, $
      color=fsc_color('firebrick',101), symsize=symsize, $
      thick=5.0
    im_legend, 'Metal-Rich HII Regions', /left, /top, box=0, $
      psym=symcat(6,thick=5.0), color='firebrick', $
      charsize=1.7

; excitation arrow    
    arrow, -0.2, 0.55, -0.2, 0.8, /data, thick=6.0, $
      hsize=-0.15, hthick=6.0
    xyouts, -0.1, 0.65, 'increasing excitation', /data, $
      orientation=90, align=0.5, charsize=1.4
; metallicity arrow    
    arrow, 1.05, 0.1, 0.7, 0.1, /data, thick=6.0, $
      hsize=-0.15, hthick=6.0
    xyouts, 0.875, 0.05, 'increasing metallicity', /data, $
      orientation=0, align=0.5, charsize=1.4
    im_plotconfig, /psclose

stop    
    
return
end
