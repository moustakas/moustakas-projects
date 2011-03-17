pro test_usno_cuts
; jm10feb08ucsd - 

    agespath = ages_path(/analysis)

;   usno = mrdfits(agespath+'catalog.usno.fits.gz',1)
;   rad = 10.0+2.5*(15.0-usno.r_usno)
;   flag = where((usno.bstar eq 1) and (usno.dstar lt rad))
    
; remove duplicates    
    ages1 = mrdfits(agespath+'catalog.cat.noguidestars.fits.gz',1)
    weight = mrdfits(agespath+'catalog.spectweight.fits.gz',1)
    ausno = mrdfits(agespath+'catalog.usno.fits.gz',1)

    usno = mrdfits(agespath+'ages_usno.fits.gz',1)
    ing = spheregroup(usno.ra,usno.dec,30.0/3600.0,$
      firstg=firstg,multg=multg,nextg=nextg)
    firstg = firstg[0L:max(ing)]
    usno = usno[firstg]
    spherematch, usno.ra, usno.dec, 15.0D*ages1.ra, $
      ages1.dec, 3.0/3600.0, m1, m2, d12
    usno1 = im_empty_structure(usno[0],empty_value=-999.0,ncopies=n_elements(ages1))
    usno1[m2] = usno[m1]
    
    win1 = ages_isin_window(15.0D*ages1.ra,ages1.dec)
    ww = where((win1 eq 0) and (weight.main_weight eq 1) and $
      (weight.spec_yesno ge 1) and (weight.z_yesno eq 1))
    struct_print, usno1[ww]
    struct_print, ausno[ww]
    

stop
; --------------------------------------------------    

    ages1 = mrdfits(agespath+'catalog.cat.noguidestars.fits.gz',1)
    usno1 = mrdfits(agespath+'catalog.usno.fits.gz',1)
    spherematch, usno.ra, usno.dec, 15.0D*ages1.ra, $
      ages1.dec, 3.0/3600.0, m1, m2, d12
    ww = where(abs(usno[m1].rmag-usno1[m2].r_usno) gt 5.0)

stop    
    
    win1 = ages_isin_window(15.0D*ages1[m2].ra,ages1[m2].dec)
    
    
stop    
    
; --------------------------------------------------
    
    agespath = ages_path(/analysis)

    allages = mrdfits(agespath+'catalog.cat.noguidestars.fits.gz',1)
    allweight = mrdfits(agespath+'catalog.spectweight.fits.gz',1)
    allusno = mrdfits(agespath+'catalog.usno.fits.gz',1)
    allcodes = mrdfits(agespath+'catalog.codes.fits.gz',1)

    allweight = mrdfits(agespath+'catalog.spectweight.fits.gz',1)
    keep = where((allweight.main_weight eq 1) and (allweight.spec_yesno ge 1) and $
      (allweight.z_yesno eq 1))
    weight = allweight[keep]
    ages = mrdfits(agespath+'catalog.cat.noguidestars.fits.gz',1,rows=keep)
    codes = mrdfits(agespath+'catalog.codes.fits.gz',1,rows=keep)
    usno = mrdfits(agespath+'catalog.usno.fits.gz',1,rows=keep)

;   win1 = ages_isin_window(15.0D*ages.ra,ages.dec)
;   check = where(win1 eq 0)
    
; test the bright-star selection
    maxis = im_array(5.0,25.0,0.02)
    im_plotconfig, 0, pos, psfile='test_usno_cuts.ps'
    djs_plot, [0], [0], /nodata, /ylog, xrange=[7,20], $
      yrange=[2,200], ps=6, sym=0.5, xsty=1, ysty=1, $
      xtitle='R_{USNO} (mag)', ytitle='DSTAR (distance to nearest bright star, arcsec)', $
      position=pos
    djs_oplot, allusno.mstar, allusno.dstar, ps=6, sym=0.2, color='orange'
    bstar = where(codes.bstar,comp=ok)
    djs_oplot, usno[bstar].mstar, usno[bstar].dstar, ps=6, sym=0.5
    djs_oplot, maxis, 20.0+5.0*(15.0-maxis), color='red'
    djs_oplot, maxis, 10.0+2.5*(15.0-maxis), color='blue'
    these = where(usno[bstar].dstar gt (20.0+5.0*(15.0-usno[bstar].mstar)))
    djs_oplot, usno[bstar[these]].mstar, usno[bstar[these]].dstar, $
      ps=6, sym=0.5, color='blue'
    legend, ['main_weight=1','bstar=1'], /left, /top, box=0, margin=0
    im_legend, ['20"+5"(15-R_{USNO})','10"+2.5"(15-R_{USNO})'], /left, /bottom, box=0, $
      color=['red','blue'], line=[0,0], pspacing=1.4, margin=0
    im_plotconfig, /psclose

    djs_oplot, allusno[ww].mstar, allusno[ww].dstar, $
      ps=6, sym=0.3, color='green'

    ww = where(allcodes.galaxy eq 0)
    bstar = where((allusno[ww].bstar eq 1),comp=ok)
    djs_plot, [0], [0], /nodata, xr=[5,20], yr=[0.1,200], $
      xsty=1, ysty=1, /ylog
    djs_oplot, allusno[ww[bstar]].mstar, allusno[ww[bstar]].dstar, $
      ps=6, sym=0.3, color='orange'
    djs_oplot, allusno[ww[ok]].mstar, allusno[ww[ok]].dstar, $
      ps=6, sym=0.3, color='purple'

    djs_oplot, maxis, 20.0+5.0*(15.0-maxis), color='red'
    djs_oplot, maxis, 10.0+2.5*(15.0-maxis), color='blue'
    djs_oplot, [15,15], 10^!y.crange, line=0


    bigrad = 20.0+5.0*(15.0-allusno.mstar)

    rad = 10.0+2.5*(15.0-allusno.mstar)
    these = where((allusno.dstar lt rad) and (allusno.mstar lt 90.0) and $
      (allusno.mstar gt 0.0) and (allusno.bstar eq 1),nstar)
    radius = rad[these]
    ra = 15.0D*allages[these].ra+allusno[these].dra/3600.0
    dec = allages[these].dec+allusno[these].ddec/3600.0
    djs_plot, 15.0D*allages.ra, allages.dec, ps=3, xsty=3, ysty=3
    for ii=0L,nstar-1 do tvcircle, 10*radius[ii]/3600.0, ra[ii], dec[ii], $
      color=djs_icolor('yellow'), /data

stop    
    
; --------------------------------------------------
    
    
    win1 = ages_isin_window(15.0D*ages.ra,ages.dec)
    struct_print, usno[where(win1 eq 0)]

;   cut1 = where((codes.bstar eq 1) and (codes.galaxy eq 1))
;   struct_print, codes[cut1]
;   win2 = ages_isin_window(15.0D*ages[cut1].ra,ages[cut1].dec)
    
    djs_plot, [0], [0], /nodata, /ylog, xrange=[7,20], $
      yrange=[0.1,200], ps=6, sym=0.3, xsty=3, ysty=3, $
      'R_{USNO} (mag)', ytitle='Distance to Bright Star (arcsec)'

    bstar = where(codes.bstar,comp=ok)
    djs_oplot, usno[bstar].r_usno, usno[bstar].dstar, ps=6, sym=0.3
    djs_oplot, usno[ok].r_usno, usno[ok].dstar, ps=6, sym=0.3, color='orange'
    djs_oplot, usno[where(win1 eq 1)].r_usno, usno[where(win1 eq 1)].dstar, $
      ps=6, sym=0.3, color='purple'
    djs_oplot, usno[where(win1 eq 0)].r_usno, usno[where(win1 eq 0)].dstar, $
      ps=6, sym=0.6, color='green'
    djs_oplot, maxis, 20.0+5.0*(15.0-maxis), color='red'
    djs_oplot, maxis, 10.0+2.5*(15.0-maxis), color='blue'

stop    

    cut1 = where(usno.bstar,ncut1)
    djs_plot, usno[cut1].r_usno, usno[cut1].dstar, $
      ps=6, sym=0.3, xsty=3, ysty=3, /ylog, xrange=[0,20]
    djs_oplot, maxis, 20.0+5.0*(15.0-maxis), color='red'

    these = where(usno[cut1].dstar gt (20.0+5.0*(15.0-usno[cut1].r_usno)))
    djs_oplot, usno[cut1[these]].r_usno, usno[cut1[these]].dstar, $
      ps=6, sym=0.3, color='blue'

    struct_print, codes[cut1[these]]
    struct_print, weight[cut1[these]]
    

stop    
    
return
end
