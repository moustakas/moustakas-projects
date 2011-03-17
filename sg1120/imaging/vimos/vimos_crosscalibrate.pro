######## OBSOLETE ########

pro vimos_crosscalibrate
; jm09mar24nyu - cross-calibrate the VIMOS photometry using the LDSS3
;   gr photometry

    common allcat, ldss3_gcat1, ldss3_rcat1, vimos_bcat1, vimos_vcat1, vimos_rcat1

    ldss3_catpath = ldss3_path(/cat)
    vimos_catpath = vimos_path(/cat)
    ldss3_version = ldss3_catalogs_version()
    vimos_version = vimos_catalogs_version()

    if (n_elements(ldss3_gcat1) eq 0L) then begin
       ldss3_gcat1 = rsex(ldss3_catpath+'sg1120_g_'+ldss3_version+'.chi2.cat')
       ldss3_rcat1 = rsex(ldss3_catpath+'sg1120_r_'+ldss3_version+'.chi2.cat')
       vimos_bcat1 = rsex(vimos_catpath+'sg1120_B_'+vimos_version+'.chi2.cat')
       vimos_vcat1 = rsex(vimos_catpath+'sg1120_V_'+vimos_version+'.chi2.cat')
       vimos_rcat1 = rsex(vimos_catpath+'sg1120_R_'+vimos_version+'.chi2.cat')
    endif

; read the SDSS catalog
    sdss1 = mrdfits(ldss3_path(/mosaics)+'sdss.cas.dr7.fits.gz',1)
    keep = where((sdss1.type eq 6) and (sdss1.gerr lt 0.1) and $
      (sdss1.rerr lt 0.1) and (sdss1.r gt 18.5) and (sdss1.r lt 23.0))
    sdss = sdss1[keep]

; crossmatch
    spherematch, ldss3_gcat1.xwin_world, ldss3_gcat1.ywin_world, $
      sdss.ra, sdss.dec, 2.0/3600.0, m1, m2
    ldss3_gcat = ldss3_gcat1[m1]
    ldss3_rcat = ldss3_rcat1[m1]

    spherematch, ldss3_gcat.xwin_world, ldss3_gcat.ywin_world, $
      vimos_vcat1.xwin_world, vimos_vcat1.ywin_world, 2.0/3600.0, $
      mm1, mm2

    good = where((ldss3_gcat[mm1].flags eq 0) and (ldss3_rcat[mm1].flags eq 0) and $
      (vimos_bcat1[mm2].flags eq 0) and (vimos_vcat1[mm2].flags eq 0) and $
      (vimos_rcat1[mm2].flags eq 0))

; correct for Galactic extinction (not negligible!)    
    euler, ldss3_gcat[mm1[good]].xwin_world, $
      ldss3_gcat[mm1[good]].ywin_world, gl, gb, 1
    ebv = dust_getval(gl,gb,/interp)
    gdust = ebv*k_lambda(k_lambda_eff(filterlist='sdss_g0.par'),/odonnell)
    rdust = ebv*k_lambda(k_lambda_eff(filterlist='sdss_r0.par'),/odonnell)
    bdust = ebv*k_lambda(k_lambda_eff(filterlist='vimos_B.par'),/odonnell)
    vdust = ebv*k_lambda(k_lambda_eff(filterlist='vimos_V.par'),/odonnell)
    rdust = ebv*k_lambda(k_lambda_eff(filterlist='vimos_R.par'),/odonnell)

;need to correct for reddening!!    
    
    g = ldss3_gcat[mm1[good]].mag_auto-gdust
    r = ldss3_rcat[mm1[good]].mag_auto-rdust
    bigb = vimos_bcat1[mm2[good]].mag_auto-bdust+k_vega2ab(filterlist='vimos_B.par',/kurucz,/silent)
    bigv = vimos_vcat1[mm2[good]].mag_auto-vdust+k_vega2ab(filterlist='vimos_V.par',/kurucz,/silent)
    bigr = vimos_rcat1[mm2[good]].mag_auto-rdust+k_vega2ab(filterlist='vimos_R.par',/kurucz,/silent)

; color transformations from Lupton+05; see
; http://www.sdss.org/dr7/algorithms/sdssUBVRITransform.html 
    
;    B = u - 0.8116*(u - g) + 0.1313;  sigma = 0.0095
;    B = g + 0.3130*(g - r) + 0.2271;  sigma = 0.0107
;  
;    V = g - 0.2906*(u - g) + 0.0885;  sigma = 0.0129
;    V = g - 0.5784*(g - r) - 0.0038;  sigma = 0.0054
;  
;    R = r - 0.1837*(g - r) - 0.0971;  sigma = 0.0106
;    R = r - 0.2936*(r - i) - 0.1439;  sigma = 0.0072
;  
;    I = r - 1.2444*(r - i) - 0.3820;  sigma = 0.0078
;    I = i - 0.3780*(i - z)  -0.3974;  sigma = 0.006

    bnew = g+0.3130*(g-r)+0.2271
    vnew = g-0.5784*(g-r)-0.0038
    rnew = r-0.1837*(g-r)-0.0971

    sigrej = 2.5
    bstats = im_stats(bnew-bigb,sigrej=sigrej)
    vstats = im_stats(vnew-bigv,sigrej=sigrej)
    rstats = im_stats(rnew-bigr,sigrej=sigrej)

    im_plotconfig, psfile=vimos_catpath+'vimos_crosscalibrate.ps'
    plotsym, 0, 1.6, /fill
    xrange = [17.5,22]
    yrange = [-1,1]
; B-band
    djs_plot, r, bnew-bigb, psym=8, xrange=xrange, yrange=yrange, $
      xsty=1, ysty=1, xtitle='r_{AB} [LDSS3]', ytitle='B_{LDSS3}-B_{VIMOS}', $
      charsize=2.0
    djs_oplot, !x.crange, bstats.mean_rej*[1,1], line=0, thick=3
    djs_oplot, !x.crange, (bstats.mean_rej+bstats.sigma_rej)*[1,1], line=5, thick=3
    djs_oplot, !x.crange, (bstats.mean_rej-bstats.sigma_rej)*[1,1], line=5, thick=3
    im_legend, '\Delta_{B} = '+im_string_stats(bnew-bigb,type=1,sigrej=sigrej), $
      /right, /top, box=0, charsize=1.8
; V-band
    djs_plot, r, vnew-bigv, psym=8, xrange=xrange, yrange=yrange, $
      xsty=1, ysty=1, xtitle='r_{AB} [LDSS3]', ytitle='V_{LDSS3}-V_{VIMOS}', $
      charsize=2.0
    djs_oplot, !x.crange, vstats.mean_rej*[1,1], line=0, thick=3
    djs_oplot, !x.crange, (vstats.mean_rej+vstats.sigma_rej)*[1,1], line=5, thick=3
    djs_oplot, !x.crange, (vstats.mean_rej-vstats.sigma_rej)*[1,1], line=5, thick=3
    im_legend, '\Delta_{V} = '+im_string_stats(vnew-bigv,type=1,sigrej=sigrej), $
      /right, /top, box=0, charsize=1.8
; R-band
    djs_plot, r, rnew-bigr, psym=8, xrange=xrange, yrange=yrange, $
      xsty=1, ysty=1, xtitle='r_{AB} [LDSS3]', ytitle='R_{LDSS3}-R_{VIMOS}', $
      charsize=2.0
    djs_oplot, !x.crange, rstats.mean_rej*[1,1], line=0, thick=3
    djs_oplot, !x.crange, (rstats.mean_rej+rstats.sigma_rej)*[1,1], line=5, thick=3
    djs_oplot, !x.crange, (rstats.mean_rej-rstats.sigma_rej)*[1,1], line=5, thick=3
    im_legend, '\Delta_{R} = '+im_string_stats(rnew-bigr,type=1,sigrej=sigrej), $
      /right, /top, box=0, charsize=1.8
    im_plotconfig, /psclose

; final zeropoints:
    bzero = bstats.median_rej
    vzero = vstats.median_rej
    rzero = rstats.median_rej
    splog, '##################################################'
    splog, 'Final zeropoints (BVR):'
    print, bzero, vzero, rzero
    
; plot the stellar locus and compare against Pickles+98    

    p98 = read_98pickles()
    filterlist = 'vimos_'+['B','V','R']+'.par'
    vega2ab = k_vega2ab(filterlist=filterlist,/kurucz,/silent)
    lambda = k_lambda_to_edges(p98[0].wave)
    maggies = fltarr(n_elements(filterlist),n_elements(p98))
    for ii = 0L, n_elements(p98)-1L do maggies[*,ii] = $
      k_project_filters(lambda,p98[ii].flux,filterlist=filterlist,/silent)
    bv_p98 = -2.5*alog10(maggies[0,*]/maggies[1,*]) - (vega2ab[0]-vega2ab[1]) ; Vega
    vr_p98 = -2.5*alog10(maggies[1,*]/maggies[2,*]) - (vega2ab[1]-vega2ab[2]) ; Vega

    im_plotconfig, psfile=vimos_catpath+'vimos_BVR_stellar_locus_calibrated.ps'
    djs_plot, [0], [0], /nodata, /xsty, /ysty, xrange=[-0.3,1.6], $
      yrange=[-0.6,2.3], xtitle='V - R', ytitle='B - V', charsize=2.0
    djs_oplot, (bigv+vzero)-(bigr-rzero), (bigb+bzero)-(bigv-vzero), psym=6, $
      color='red', sym=2.0, thick=4.0
    djs_oplot, vr_p98, bv_p98, psym=6, thick=4.0
    im_legend, ['VIMOS stars','Pickles+98'], /left, /top, box=0, $
      psym=6, color=['red',''], charsize=1.8
    im_plotconfig, /psclose
    
    
stop    
    
return
end
    
