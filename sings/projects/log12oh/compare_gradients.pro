pro compare_gradients
; jm08jun27nyu - compare my gradients with Zaritsky et al. and
; Pilyugin+04 

    metpath = sings_path(/projects)+'log12oh/'
    sings = mrdfits(metpath+'sings_log12oh_hiiregions_v8.0.fits.gz',1)
    sings = sings[where(sings.gradient_flag)]

    z94 = rsex(hiiregions_path()+'94zaritsky/spiral_gradients.dat')
    p04 = rsex(getenv('CATALOGS_DIR')+'/04pilyugin/04pilyugin_table1.dat')

    psfile = metpath+'compare_gradients.ps'
    im_plotconfig, 0, pos, psfile=psfile, xmargin=[1.3,0.2]

; Zaritsky+94    
    djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
      xtitle='Gradient (dex \rho_{25}^{-1}) [Moustakas+10, KK04]', $
      ytitle='Gradient (dex \rho_{25}^{-1}) [Zaritsky+94]', $
      xrange=[-1.5,0.3], yrange=[-1.5,0.3]
    djs_oplot, !x.crange, !y.crange, line=0
    
    match, strtrim(sings.galaxy,2), strtrim(z94.galaxy,2), m1, m2
    niceprint, sings[m1].galaxy, z94[m2].galaxy
    oploterror, sings[m1].hii_kk04_slope[0], z94[m2].gradient, $
      sings[m1].hii_kk04_slope[1], z94[m2].gradient_err, $
      errthick=0.1, psym=symcat(16), symsize=1.5

; Pilyugin+04
    djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
      xtitle='Gradient (dex \rho_{25}^{-1}) [Moustakas+10, PT05]', $
      ytitle='Gradient (dex \rho_{25}^{-1}) [Pilyugin+04]', $
      xrange=[-1.5,0.3], yrange=[-1.5,0.3]
    djs_oplot, !x.crange, !y.crange, line=0
    
    match, strtrim(sings.galaxy,2), strtrim(p04.galaxy,2), m1, m2
    niceprint, sings[m1].galaxy, p04[m2].galaxy
    oploterror, sings[m1].hii_kk04_slope[0], p04[m2].ohgradient_r25, $
      sings[m1].hii_kk04_slope[1], p04[m2].ohgradient_r25*0.0, $
      errthick=0.1, psym=symcat(16), symsize=1.5
    for ii = 0, n_elements(m1)-1 do xyouts, sings[m1[ii]].hii_kk04_slope[0], $
      p04[m2[ii]].ohgradient_r25+0.05, strtrim(sings[m1[ii]].galaxy,2), align=0.0, $
      charsize=1.1

    im_plotconfig, psfile=psfile, /psclose, /gzip

return
end
    
