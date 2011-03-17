


; ------------------------------------------------------------
; SFR(Ha) vs SFR(Hb) & SFR([O II])
; ------------------------------------------------------------

; need to apply aperture corrections!!    

    psname = 'ages_sfr_ha_vs_sfr_hb_sfr_oii'
    im_openclose, pspath+psname, postscript=postscript, encapsulated=encapsulated, xsize=6.5, ysize=8.5

    pagemaker, nx=1, ny=2, height=3.5*[1.0,1.0], width=5.0, xmargin=[1.1,0.4], $
      ymargin=[0.4,1.1], xspace=0, yspace=0, xpage=6.5, ypage=8.5, $
      position=pos, /normal

    xrange = sfrharange
    yrange = xrange

    xtitle = 'log \psi(H\alpha) ['+sfr_units()+']'
    ytitle1 = 'log \psi(H\beta) ['+sfr_units()+']'
    ytitle2 = 'log \psi([O II]) ['+sfr_units()+']'    

    indx = where((agesancillary.m_b gt -900.0) and (agesancillary.kcorr_mass gt -900.0) and $
      (agesdust.zstrong_ew_12oh_kk04 gt -900) and (strmatch(agesdust.r23branch_ew_kk04,'*U*') eq 1B) and $
      (agesancillary.fluxing_problem eq 0B) and (agesnodust.sfr_h_alpha gt -900.0) and $
      (agesdust.h_beta_lum[0] gt -900.0) and (agesdust.oii_3727_lum[0] gt -900.0),nindx)
    sf = where(strtrim(agesdust[indx].final_class,2) eq 'BPT_HII',nsf)
    unk = where(strtrim(agesdust[indx].final_class,2) eq '?',nunk)
    
    hasfr = agesnodust[indx].sfr_h_alpha
    hasfrerr = agesnodust[indx].sfr_h_alpha_err

; SFR(Ha) vs SFR(Hb)

    lhbobs = agesdust[indx].h_beta_lum[0] + alog10(lsun)
    loglb = agesancillary[indx].b_lum
    hbsfr = hb_sfr(loglb,lhbobs,sfr_err=hbsfrerr,/log)

    stats = im_stats(hbsfr-hasfr,/verbose,/baremin)
    xstr = '\Delta'+'log(\psi)='+strtrim(string(stats.median_rej,format='(F12.2)'),2)+'\pm'+$
      strtrim(string(stats.sigma_rej,format='(F12.2)'),2)

    ages_lineplot, hasfr, hbsfr, hasfrerr, hbsfrerr, plottype=1, postscript=postscript, /bin2d, $
      xtitle='', ytitle=ytitle1, xrange=xrange, yrange=yrange, $
      charsize=charsize_6, xtickname=replicate(' ',10), position=pos[*,0], yminor=5
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick

    legend, '(a)', /left, /top, box=0, charsize=charsize_6, charthick=postthick
    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_3, charthick=postthick

; SFR(Ha) vs SFR([OII])
    
    loiiobs = agesdust[indx].oii_3727_lum[0] + alog10(lsun)
    loglb = agesancillary[indx].b_lum
    oiisfr = oii_sfr(loglb,loiiobs,sfr_err=oiisfrerr,/log)

    stats = im_stats(oiisfr-hasfr,/verbose,/baremin,/no_head)
    xstr = '\Delta'+'log(\psi)='+strtrim(string(stats.median_rej,format='(F12.2)'),2)+'\pm'+$
      strtrim(string(stats.sigma_rej,format='(F12.2)'),2)
    
    ages_lineplot, hasfr, oiisfr, hasfrerr, oiisfrerr, plottype=1, postscript=postscript, $
      /bin2d, /noerase, xtitle=xtitle, ytitle=ytitle2, xrange=xrange, yrange=yrange, $
      charsize=charsize_6, position=pos[*,1], yminor=5
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick

    legend, '(b)', /left, /top, box=0, charsize=charsize_6, charthick=postthick
    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_3, charthick=postthick

    im_openclose, postscript=postscript, /close    


; ------------------------------------------------------------
; 4-panel SFR vs M_B, in bins of redshift
; ------------------------------------------------------------

    psname = 'ages_mb_vs_sfr_4panel_redshift'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=6.5, encapsulated=encapsulated

    pagemaker, nx=2, ny=2, yspace=0, xspace=0, width=3.5*[1,1], height=2.5*[1,1], $
      xmargin=[1.1,0.4], ymargin=[0.4,1.1], xpage=8.5, ypage=6.5, position=pos, /normal

    indx = where((agesancillary.m_b gt -900.0) and (agesancillary.kcorr_mass gt -900.0) and $
      (agesdust.zstrong_ew_12oh_kk04 gt -900) and (strmatch(agesdust.r23branch_ew_kk04,'*U*') eq 1B) and $
      (agesancillary.fluxing_problem eq 0B),nindx)

    mb = agesancillary[indx].m_b
    mberr = agesancillary[indx].m_b_err
    zobj = agesdust[indx].z_obj
    oh = agesdust[indx].zstrong_ew_12oh_kk04
    oherr = agesdust[indx].zstrong_ew_12oh_kk04_err

    lhbobs = agesdust[indx].h_beta_lum[0] + alog10(lsun)
    loglb = agesancillary[indx].b_lum
    sfr = hb_sfr(loglb,lhbobs,sfr_err=sfrerr,/log)

    intindx = where((intdust.m_b_obs gt -900.0) and (intdust.zstrong_ew_12oh_kk04 gt -900.0) and $
      (strmatch(intdust.r23branch_ew_kk04,'*U*') eq 1B),nintindx)

    mbint = intdust[intindx].m_b_obs
    mbinterr = intdust[intindx].m_b_obs_err
    
    lhbobsint = intdust[intindx].h_beta_lum[0] + alog10(lsun)
    loglbint = intdust[intindx].b_lum_obs
    sfrint = hb_sfr(loglbint,lhbobsint,sfr_err=sfrinterr,/log)

    xtitle = 'M_{B} [mag]'
    ytitle = 'log \psi(H\beta) ['+sfr_units()+']'

    xrange = mbrange7
    yrange = sfrrange2

    yerrpos = 0.6
    
; ##########
; z<0.2
; ##########

    bin1 = where((zobj le 0.2),nbin1)
    hi = where(oh[bin1] gt agesohcut,nhi,comp=lo,ncomp=nlo)

    x = mb[bin1]
    xerr = mberr[bin1]
    y = sfr[bin1]
    yerr = sfrerr[bin1]

    medyerr = median(yerr)
    medxerr = median(xerr)

; draw "lo" before "hi" because of the colors and the point density    
    
    ages_lineplot, x[lo], y[lo], xerr[lo], yerr[lo], plottype=1, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, $
      position=pos[*,0], charsize=charsize_6, xtickname=replicate(' ',10), $
      agespsize=agesohlopsize, agessym=agesohlosym, agescolor=agesohlocolor
    if (nhi ne 0L) then ages_lineplot, x[hi], y[hi], xerr[hi], yerr[hi], $
      /overplot, agespsize=agesohhipsize, agessym=agesohhisym, agescolor=agesohhicolor
    djs_oplot, -18*[1,1], !y.crange, line=1, thick=postthick2
;   oploterror, !x.crange[0]+yerrpos, !y.crange[0]+yerrpos, medxerr, $
;     medyerr, ps=3, /data, thick=postthick3, errthick=postthick
    legend, '0.0<z<0.2', /left, /top, box=0, charsize=charsize_2, charthick=postthick

;   atlas1d_lineplot, mbint, sfrint, mbinterr, sfrinterr, /overplot, $
;     postscript=postscript, atlascolor=intcolor2, atlaspsize=intpsize, atlassym=intsym

; ##########
; 0.2<z<0.4
; ##########

    bin1 = where((zobj gt 0.2) and (zobj le 0.4),nbin1)
    hi = where(oh[bin1] gt agesohcut,nhi,comp=lo,ncomp=nlo)

    x = mb[bin1]
    xerr = mberr[bin1]
    y = sfr[bin1]
    yerr = sfrerr[bin1]

    medyerr = median(yerr)
    medxerr = median(xerr)

; draw "lo" before "hi" because of the colors and the point density    
    
    ages_lineplot, x[lo], y[lo], xerr[lo], yerr[lo], plottype=1, postscript=postscript, $
      xtitle='', ytitle='', xrange=xrange, yrange=yrange, $
      position=pos[*,1], charsize=charsize_6, xtickname=replicate(' ',10), $
      ytickname=replicate(' ',10), /noerase, $
      agespsize=agesohlopsize, agessym=agesohlosym, agescolor=agesohlocolor
    if (nhi ne 0L) then ages_lineplot, x[hi], y[hi], xerr[hi], yerr[hi], $
      /overplot, agespsize=agesohhipsize, agessym=agesohhisym, agescolor=agesohhicolor
    djs_oplot, -19.0*[1,1], !y.crange, line=1, thick=postthick2
    legend, '0.2<z<0.4', /left, /top, box=0, charsize=charsize_2, charthick=postthick
;   oploterror, !x.crange[0]+yerrpos, !y.crange[0]+yerrpos, medxerr, $
;     medyerr, ps=3, /data, thick=postthick3, errthick=postthick

;   atlas1d_lineplot, mbint, sfrint, mbinterr, sfrinterr, /overplot, $
;     postscript=postscript, atlascolor=intcolor2, atlaspsize=intpsize, atlassym=intsym

; ##########
; 0.4<z<0.6
; ##########

    bin1 = where((zobj gt 0.4) and (zobj le 0.6),nbin1)
    hi = where(oh[bin1] gt agesohcut,nhi,comp=lo,ncomp=nlo)

    x = mb[bin1]
    xerr = mberr[bin1]
    y = sfr[bin1]
    yerr = sfrerr[bin1]

    medyerr = median(yerr)
    medxerr = median(xerr)

    ages_lineplot, x[hi], y[hi], xerr[hi], yerr[hi], plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, $
      position=pos[*,2], charsize=charsize_6, /noerase, $
      agespsize=agesohhipsize*1.4, agessym=agesohhisym, agescolor=agesohhicolor
    if (nlo ne 0L) then ages_lineplot, x[lo], y[lo], xerr[lo], yerr[lo], $
      /overplot, agespsize=agesohlopsize*1.2, agessym=agesohlosym, agescolor=agesohlocolor
    djs_oplot, -20.5*[1,1], !y.crange, line=1, thick=postthick2
    legend, '0.4<z<0.6', /left, /top, box=0, charsize=charsize_2, charthick=postthick
;   oploterror, !x.crange[0]+yerrpos, !y.crange[0]+yerrpos, medxerr, $
;     medyerr, ps=3, /data, thick=postthick3, errthick=postthick

;   atlas1d_lineplot, mbint, sfrint, mbinterr, sfrinterr, /overplot, $
;     postscript=postscript, atlascolor=intcolor2, atlaspsize=intpsize, atlassym=intsym

; ##########
; 0.6<z<0.8
; ##########

    bin1 = where((zobj gt 0.6) and (zobj le max(zobj)),nbin1)
    hi = where(oh[bin1] gt agesohcut,nhi,comp=lo,ncomp=nlo)

    x = mb[bin1]
    xerr = mberr[bin1]
    y = sfr[bin1]
    yerr = sfrerr[bin1]

    medyerr = median(yerr)
    medxerr = median(xerr)

    ages_lineplot, x[hi], y[hi], xerr[hi], yerr[hi], plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle='', xrange=xrange, yrange=yrange, $
      position=pos[*,3], charsize=charsize_6, /noerase, ytickname=replicate(' ',10), $
      agespsize=agesohhipsize*1.4, agessym=agesohhisym, agescolor=agesohhicolor
    if (nlo ne 0L) then ages_lineplot, x[lo], y[lo], xerr[lo], yerr[lo], $
      /overplot, agespsize=agesohlopsize*1.2, agessym=agesohlosym, agescolor=agesohlocolor
    djs_oplot, -21.5*[1,1], !y.crange, line=1, thick=postthick2
    legend, '0.6<z<0.8', /left, /top, box=0, charsize=charsize_2, charthick=postthick
;   oploterror, !x.crange[0]+yerrpos, !y.crange[0]+yerrpos, medxerr, $
;     medyerr, ps=3, /data, thick=postthick3, errthick=postthick

;   atlas1d_lineplot, mbint, sfrint, mbinterr, sfrinterr, /overplot, $
;     postscript=postscript, atlascolor=intcolor2, atlaspsize=intpsize, atlassym=intsym

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; 4-panel SSFR vs mass, in bins of redshift
; ------------------------------------------------------------

    psname = 'ages_mass_vs_ssfr_4panel_redshift'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=6.5, encapsulated=encapsulated

    pagemaker, nx=2, ny=2, yspace=0, xspace=0, width=3.5*[1,1], height=2.5*[1,1], $
      xmargin=[1.1,0.4], ymargin=[0.4,1.1], xpage=8.5, ypage=6.5, position=pos, /normal

    indx = where((agesancillary.m_b gt -900.0) and (agesancillary.kcorr_mass gt -900.0) and $
      (agesdust.zstrong_ew_12oh_kk04 gt -900) and (strmatch(agesdust.r23branch_ew_kk04,'*U*') eq 1B) and $
      (agesancillary.fluxing_problem eq 0B),nindx)

    mass = agesancillary[indx].kcorr_mass
    masserr = mass*0.0+masssyserr
    zobj = agesdust[indx].z_obj
    oh = agesdust[indx].zstrong_ew_12oh_kk04
    oherr = agesdust[indx].zstrong_ew_12oh_kk04_err

    lhbobs = agesdust[indx].h_beta_lum[0] + alog10(lsun)
    loglb = agesancillary[indx].b_lum
    sfr = hb_sfr(loglb,lhbobs,sfr_err=sfrerr,/log)
    ssfr = sfr - mass + alog10(1E9) ; [Gyr^-1]
    ssfrerr = sqrt(sfrerr^2.0 + masserr^2.0)

    intindx = where((intdust.m_b_obs gt -900.0) and (intdust.kcorr_mass gt -900.0) and $
      (intdust.zstrong_ew_12oh_kk04 gt -900.0) and (strmatch(intdust.r23branch_ew_kk04,'*U*') eq 1B),nintindx)

    massint = intdust[intindx].kcorr_mass
    massinterr = massint*0.0+masssyserr
    
    lhbobsint = intdust[intindx].h_beta_lum[0] + alog10(lsun)
    loglbint = intdust[intindx].b_lum_obs
    sfrint = hb_sfr(loglbint,lhbobsint,sfr_err=sfrerrint,/log)
    ssfrint = sfrint - massint + alog10(1E9) ; [Gyr^-1]
    ssfrinterr = sqrt(sfrerrint^2.0 + massinterr^2.0)

    xtitle = 'log M [M'+sunsymbol()+']'
    ytitle = 'log (\psi/M) [Gyr^{-1}]'

    xrange = massrange3
    yrange = ssfrrange

    yerrpos = 0.6
    
    ssfrconst1 = alog10(0.1) - massaxis + alog10(1E9)  
    ssfrconst2 = alog10(1.0) - massaxis + alog10(1E9)  ; sfr = 1 M_sun/yr [Gyr^-1]
    ssfrconst3 = alog10(10.0) - massaxis + alog10(1E9)
    ssfrconst4 = alog10(100.0) - massaxis + alog10(1E9)

; ##########
; z<0.2
; ##########

    bin1 = where((zobj le 0.2),nbin1)
    hi = where(oh[bin1] gt agesohcut,nhi,comp=lo,ncomp=nlo)

    x = mass[bin1]
    xerr = masserr[bin1]
    y = ssfr[bin1]
    yerr = ssfrerr[bin1]

    medyerr = median(yerr)
    medxerr = median(xerr)

; draw "lo" before "hi" because of the colors and the point density    
    
    ages_lineplot, x[lo], y[lo], xerr[lo], yerr[lo], plottype=1, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, $
      position=pos[*,0], charsize=charsize_6, xtickname=replicate(' ',10), $
      agespsize=agesohlopsize, agessym=agesohlosym, agescolor=agesohlocolor
    if (nhi ne 0L) then ages_lineplot, x[hi], y[hi], xerr[hi], yerr[hi], $
      /overplot, agespsize=agesohhipsize, agessym=agesohhisym, agescolor=agesohhicolor
    legend, '0.0<z<0.2', /right, /top, box=0, charsize=charsize_2, charthick=postthick
    oploterror, !x.crange[0]+yerrpos, !y.crange[0]+yerrpos, medxerr, $
      medyerr, ps=3, /data, thick=postthick3, errthick=postthick

;   atlas1d_lineplot, massint, ssfrint, massinterr, ssfrinterr, /overplot, $
;     postscript=postscript, atlascolor=intcolor2, atlaspsize=intpsize, atlassym=intsym

;   djs_oplot, massaxis, ssfrconst1, line=0, thick=postthick2
    djs_oplot, massaxis, ssfrconst2, line=1, thick=postthick2
    djs_oplot, massaxis, ssfrconst3, line=2, thick=postthick2
    djs_oplot, massaxis, ssfrconst4, line=3, thick=postthick2
;   djs_oplot, 9.2*[1,1], !y.crange, line=1, thick=postthick2

; ##########
; 0.2<z<0.4
; ##########

    bin1 = where((zobj gt 0.2) and (zobj le 0.4),nbin1)
    hi = where(oh[bin1] gt agesohcut,nhi,comp=lo,ncomp=nlo)

    x = mass[bin1]
    xerr = masserr[bin1]
    y = ssfr[bin1]
    yerr = ssfrerr[bin1]

    medyerr = median(yerr)
    medxerr = median(xerr)

; draw "lo" before "hi" because of the colors and the point density    
    
    ages_lineplot, x[lo], y[lo], xerr[lo], yerr[lo], plottype=1, postscript=postscript, $
      xtitle='', ytitle='', xrange=xrange, yrange=yrange, $
      position=pos[*,1], charsize=charsize_6, xtickname=replicate(' ',10), $
      ytickname=replicate(' ',10), /noerase, $
      agespsize=agesohlopsize, agessym=agesohlosym, agescolor=agesohlocolor
    if (nhi ne 0L) then ages_lineplot, x[hi], y[hi], xerr[hi], yerr[hi], $
      /overplot, agespsize=agesohhipsize, agessym=agesohhisym, agescolor=agesohhicolor
    legend, '0.2<z<0.4', /right, /top, box=0, charsize=charsize_2, charthick=postthick
    oploterror, !x.crange[0]+yerrpos, !y.crange[0]+yerrpos, medxerr, $
      medyerr, ps=3, /data, thick=postthick3, errthick=postthick

;   atlas1d_lineplot, massint, ssfrint, massinterr, ssfrinterr, /overplot, $
;     postscript=postscript, atlascolor=intcolor2, atlaspsize=intpsize, atlassym=intsym

;   djs_oplot, massaxis, ssfrconst1, line=0, thick=postthick2
    djs_oplot, massaxis, ssfrconst2, line=1, thick=postthick2
    djs_oplot, massaxis, ssfrconst3, line=2, thick=postthick2
    djs_oplot, massaxis, ssfrconst4, line=3, thick=postthick2
;   djs_oplot, 10.15*[1,1], !y.crange, line=1, thick=postthick2

; ##########
; 0.4<z<0.6
; ##########

    bin1 = where((zobj gt 0.4) and (zobj le 0.6),nbin1)
    hi = where(oh[bin1] gt agesohcut,nhi,comp=lo,ncomp=nlo)

    x = mass[bin1]
    xerr = masserr[bin1]
    y = ssfr[bin1]
    yerr = ssfrerr[bin1]

    medyerr = median(yerr)
    medxerr = median(xerr)

    ages_lineplot, x[hi], y[hi], xerr[hi], yerr[hi], plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, $
      position=pos[*,2], charsize=charsize_6, /noerase, $
      agespsize=agesohhipsize*1.4, agessym=agesohhisym, agescolor=agesohhicolor
    if (nlo ne 0L) then ages_lineplot, x[lo], y[lo], xerr[lo], yerr[lo], $
      /overplot, agespsize=agesohlopsize*1.2, agessym=agesohlosym, agescolor=agesohlocolor
    legend, '0.4<z<0.6', /right, /top, box=0, charsize=charsize_2, charthick=postthick
    oploterror, !x.crange[0]+yerrpos, !y.crange[0]+yerrpos, medxerr, $
      medyerr, ps=3, /data, thick=postthick3, errthick=postthick

;   atlas1d_lineplot, massint, ssfrint, massinterr, ssfrinterr, /overplot, $
;     postscript=postscript, atlascolor=intcolor2, atlaspsize=intpsize, atlassym=intsym

;   djs_oplot, massaxis, ssfrconst1, line=0, thick=postthick2
    djs_oplot, massaxis, ssfrconst2, line=1, thick=postthick2
    djs_oplot, massaxis, ssfrconst3, line=2, thick=postthick2
    djs_oplot, massaxis, ssfrconst4, line=3, thick=postthick2
;   djs_oplot, 10.65*[1,1], !y.crange, line=1, thick=postthick2

; ##########
; 0.6<z<0.8
; ##########

    bin1 = where((zobj gt 0.6) and (zobj le max(zobj)),nbin1)
    hi = where(oh[bin1] gt agesohcut,nhi,comp=lo,ncomp=nlo)

    x = mass[bin1]
    xerr = masserr[bin1]
    y = ssfr[bin1]
    yerr = ssfrerr[bin1]

    medyerr = median(yerr)
    medxerr = median(xerr)

    ages_lineplot, x[hi], y[hi], xerr[hi], yerr[hi], plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle='', xrange=xrange, yrange=yrange, $
      position=pos[*,3], charsize=charsize_6, /noerase, ytickname=replicate(' ',10), $
      agespsize=agesohhipsize*1.4, agessym=agesohhisym, agescolor=agesohhicolor
    if (nlo ne 0L) then ages_lineplot, x[lo], y[lo], xerr[lo], yerr[lo], $
      /overplot, agespsize=agesohlopsize*1.2, agessym=agesohlosym, agescolor=agesohlocolor
    legend, '0.6<z<0.8', /right, /top, box=0, charsize=charsize_2, charthick=postthick
    oploterror, !x.crange[0]+yerrpos, !y.crange[0]+yerrpos, medxerr, $
      medyerr, ps=3, /data, thick=postthick3, errthick=postthick

;   atlas1d_lineplot, massint, ssfrint, massinterr, ssfrinterr, /overplot, $
;     postscript=postscript, atlascolor=intcolor2, atlaspsize=intpsize, atlassym=intsym

;   djs_oplot, massaxis, ssfrconst1, line=0, thick=postthick2
    djs_oplot, massaxis, ssfrconst2, line=1, thick=postthick2
    djs_oplot, massaxis, ssfrconst3, line=2, thick=postthick2
    djs_oplot, massaxis, ssfrconst4, line=3, thick=postthick2
;   djs_oplot, 10.95*[1,1], !y.crange, line=1, thick=postthick2

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; 4-panel SSFR vs 12+log(O/H), in bins of redshift
; ------------------------------------------------------------

    psname = 'ages_12oh_vs_ssfr_4panel_redshift'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=6.5, encapsulated=encapsulated

    pagemaker, nx=2, ny=2, yspace=0, xspace=0, width=3.5*[1,1], height=2.5*[1,1], $
      xmargin=[1.1,0.4], ymargin=[0.4,1.1], xpage=8.5, ypage=6.5, position=pos, /normal

    indx = where((agesancillary.m_b gt -900.0) and (agesancillary.kcorr_mass gt -900.0) and $
      (agesdust.zstrong_ew_12oh_kk04 gt -900) and (strmatch(agesdust.r23branch_ew_kk04,'*U*') eq 1B) and $
      (agesancillary.fluxing_problem eq 0B),nindx)

    mass = agesancillary[indx].kcorr_mass
    masserr = mass*0.0+masssyserr
    zobj = agesdust[indx].z_obj
    oh = agesdust[indx].zstrong_ew_12oh_kk04
    oherr = agesdust[indx].zstrong_ew_12oh_kk04_err

    lhbobs = agesdust[indx].h_beta_lum[0] + alog10(lsun)
    loglb = agesancillary[indx].b_lum
    sfr = hb_sfr(loglb,lhbobs,sfr_err=sfrerr,/log)
    ssfr = sfr - mass + alog10(1E9) ; [Gyr^-1]
    ssfrerr = sqrt(sfrerr^2.0 + masserr^2.0)

    intindx = where((intdust.m_b_obs gt -900.0) and (intdust.kcorr_mass gt -900.0) and $
      (intdust.zstrong_ew_12oh_kk04 gt -900.0) and (strmatch(intdust.r23branch_ew_kk04,'*U*') eq 1B),nintindx)

    massint = intdust[intindx].kcorr_mass
    massinterr = massint*0.0+masssyserr
    ohint = intdust[intindx].zstrong_ew_12oh_kk04
    ohinterr = intdust[intindx].zstrong_ew_12oh_kk04_err
    
    lhbobsint = intdust[intindx].h_beta_lum[0] + alog10(lsun)
    loglbint = intdust[intindx].b_lum_obs
    sfrint = hb_sfr(loglbint,lhbobsint,sfr_err=sfrerrint,/log)
    ssfrint = sfrint - massint + alog10(1E9) ; [Gyr^-1]
    ssfrinterr = sqrt(sfrerrint^2.0 + massinterr^2.0)

    xtitle = '12 + log (O/H)_{EW}'
    ytitle = 'log (\psi/M) [Gyr^{-1}]'

    xrange = kk04ohrange4
    yrange = ssfrrange

    xerrpos = 0.25
    yerrpos = 0.50
    
; ##########
; z<0.2
; ##########

    bin1 = where((zobj le 0.2),nbin1)

    x = oh[bin1]
    xerr = oherr[bin1]
    y = ssfr[bin1]
    yerr = ssfrerr[bin1]

    medyerr = median(yerr)
    medxerr = median(xerr)

    ages_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, $
      position=pos[*,0], charsize=charsize_6, xtickname=replicate(' ',10), $
      agespsize=ageslitpsize, agessym=ageslitsym, agescolor=ageslitcolor
    legend, '(a) 0.0<z<0.2', /left, /top, box=0, charsize=charsize_2, charthick=postthick

    oploterror, !x.crange[0]+ohewsyserr*1.5, !y.crange[0]+yerrpos, ohewsyserr, $
      medyerr, ps=3, /data, thick=postthick3, errthick=postthick
    oploterror, !x.crange[0]+ohewsyserr*1.5+xerrpos, !y.crange[0]+yerrpos, medxerr, $
      medyerr, ps=3, /data, thick=postthick3, errthick=postthick

;   atlas1d_lineplot, ohint, ssfrint, ohinterr, ssfrinterr, /overplot, $
;     postscript=postscript, atlascolor=intcolor, atlaspsize=intpsize, atlassym=intsym

; ##########
; 0.2<z<0.4
; ##########

    bin1 = where((zobj gt 0.2) and (zobj le 0.4),nbin1)

    x = oh[bin1]
    xerr = oherr[bin1]
    y = ssfr[bin1]
    yerr = ssfrerr[bin1]

    medyerr = median(yerr)
    medxerr = median(xerr)

    ages_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle='', ytitle='', xrange=xrange, yrange=yrange, $
      position=pos[*,1], charsize=charsize_6, xtickname=replicate(' ',10), $
      ytickname=replicate(' ',10), /noerase, $
      agespsize=ageslitpsize, agessym=ageslitsym, agescolor=ageslitcolor
    legend, '(b) 0.2<z<0.4', /left, /top, box=0, charsize=charsize_2, charthick=postthick

    oploterror, !x.crange[0]+ohewsyserr*1.5, !y.crange[0]+yerrpos, ohewsyserr, $
      medyerr, ps=3, /data, thick=postthick3, errthick=postthick
    oploterror, !x.crange[0]+ohewsyserr*1.5+xerrpos, !y.crange[0]+yerrpos, medxerr, $
      medyerr, ps=3, /data, thick=postthick3, errthick=postthick

;   atlas1d_lineplot, ohint, ssfrint, ohinterr, ssfrinterr, /overplot, $
;     postscript=postscript, atlascolor=intcolor, atlaspsize=intpsize, atlassym=intsym

; ##########
; 0.4<z<0.6
; ##########

    bin1 = where((zobj gt 0.4) and (zobj le 0.6),nbin1)

    x = oh[bin1]
    xerr = oherr[bin1]
    y = ssfr[bin1]
    yerr = ssfrerr[bin1]

    medyerr = median(yerr)
    medxerr = median(xerr)

    ages_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, $
      position=pos[*,2], charsize=charsize_6, /noerase, $
      agespsize=ageslitpsize, agessym=ageslitsym, agescolor=ageslitcolor
    legend, '(c) 0.4<z<0.6', /left, /top, box=0, charsize=charsize_2, charthick=postthick

    oploterror, !x.crange[0]+ohewsyserr*1.5, !y.crange[0]+yerrpos, ohewsyserr, $
      medyerr, ps=3, /data, thick=postthick3, errthick=postthick
    oploterror, !x.crange[0]+ohewsyserr*1.5+xerrpos, !y.crange[0]+yerrpos, medxerr, $
      medyerr, ps=3, /data, thick=postthick3, errthick=postthick

;   atlas1d_lineplot, ohint, ssfrint, ohinterr, ssfrinterr, /overplot, $
;     postscript=postscript, atlascolor=intcolor, atlaspsize=intpsize, atlassym=intsym

; ##########
; 0.6<z<0.8
; ##########

    bin1 = where((zobj gt 0.6) and (zobj le max(zobj)),nbin1)

    x = oh[bin1]
    xerr = oherr[bin1]
    y = ssfr[bin1]
    yerr = ssfrerr[bin1]

    medyerr = median(yerr)
    medxerr = median(xerr)

    ages_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle='', xrange=xrange, yrange=yrange, $
      position=pos[*,3], charsize=charsize_6, /noerase, ytickname=replicate(' ',10), $
      agespsize=ageslitpsize, agessym=ageslitsym, agescolor=ageslitcolor
    legend, '(d) 0.6<z<0.8', /left, /top, box=0, charsize=charsize_2, charthick=postthick

    oploterror, !x.crange[0]+ohewsyserr*1.5, !y.crange[0]+yerrpos, ohewsyserr, $
      medyerr, ps=3, /data, thick=postthick3, errthick=postthick
    oploterror, !x.crange[0]+ohewsyserr*1.5+xerrpos, !y.crange[0]+yerrpos, medxerr, $
      medyerr, ps=3, /data, thick=postthick3, errthick=postthick

;   atlas1d_lineplot, ohint, ssfrint, ohinterr, ssfrinterr, /overplot, $
;     postscript=postscript, atlascolor=intcolor, atlaspsize=intpsize, atlassym=intsym

    im_openclose, postscript=postscript, /close    


; ------------------------------------------------------------
; 4-panel M/L vs mass, in bins of redshift
; ------------------------------------------------------------

    psname = 'ages_mass_vs_ml_4panel_redshift'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=6.5, encapsulated=encapsulated

    pagemaker, nx=2, ny=2, yspace=0, xspace=0, width=3.5*[1,1], height=2.5*[1,1], $
      xmargin=[1.1,0.4], ymargin=[0.4,1.1], xpage=8.5, ypage=6.5, position=pos, /normal

    indx = where((agesancillary.m_b gt -900.0) and (agesancillary.kcorr_mass gt -900.0) and $
      (agesdust.zstrong_ew_12oh_kk04 gt -900) and (agesancillary.fluxing_problem eq 0B),nindx)

    mass = agesancillary[indx].kcorr_mass
    masserr = mass*0.0+masssyserr
    zobj = agesdust[indx].z_obj
    oh = agesdust[indx].zstrong_ew_12oh_kk04
    oherr = agesdust[indx].zstrong_ew_12oh_kk04_err

    loglb = agesancillary[indx].b_lum
    loglberr = agesancillary[indx].b_lum_err

    ml = mass - loglb - 0.4*(mbolsun-mbsun) ; normalize to the total solar luminosity
    mlerr = sqrt(masserr + loglberr^2)

    intindx = where((intdust.kcorr_mass gt -900.0) and (intdust.m_b_obs gt -900.0) and $
      (intdust.zstrong_ew_12oh_kk04 gt -900.0) and (strmatch(intdust.r23branch_ew_kk04,'*U*') eq 1B),nintindx)
    massint = intdust[intindx].kcorr_mass
    massinterr = massint*0.0+masssyserr
    
    loglbint = intdust[intindx].b_lum_obs
    loglbinterr = intdust[intindx].b_lum_obs_err
    
    mlint = massint - loglb - 0.4*(mbolsun-mbsun) ; normalize to the total solar luminosity
    mlinterr = sqrt(massinterr + loglbinterr^2)
    
    xtitle = 'log M [M'+sunsymbol()+']'
    ytitle = 'log (M/L_{B})'

    xrange = massrange5
    yrange = mlrange

    yerrpos = 0.8
    
    lconst1 = massaxis - alog10(1E10) - 0.0/2.0 ; reference
    lconst2 = massaxis - alog10(1E10) - 1.0/2.5 ; 1 magnitude brighter
    lconst3 = massaxis - alog10(1E10) - 2.0/2.5 ; 2 magnitudes brighter

; ##########
; z<0.2
; ##########

    bin1 = where((zobj le 0.2),nbin1)
    hi = where(oh[bin1] gt agesohcut,nhi,comp=lo,ncomp=nlo)

    x = mass[bin1]
    xerr = masserr[bin1]
    y = ml[bin1]
    yerr = mlerr[bin1]

    medyerr = median(yerr)
    medxerr = median(xerr)

; draw "lo" before "hi" because of the colors and the point density    
    
    ages_lineplot, x[lo], y[lo], xerr[lo], yerr[lo], plottype=1, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, $
      position=pos[*,0], charsize=charsize_6, xtickname=replicate(' ',10), $
      agespsize=agesohlopsize, agessym=agesohlosym, agescolor=agesohlocolor
    if (nhi ne 0L) then ages_lineplot, x[hi], y[hi], xerr[hi], yerr[hi], $
      /overplot, agespsize=agesohhipsize, agessym=agesohhisym, agescolor=agesohhicolor
    legend, '0.0<z<0.2', /left, /top, box=0, charsize=charsize_2, charthick=postthick
;   oploterror, !x.crange[1]-yerrpos, !y.crange[0]+yerrpos*1.2, medxerr, $
;     medyerr, ps=3, /data, thick=postthick3, errthick=postthick

;   atlas1d_lineplot, massint, mlint, massinterr, mlinterr, /overplot, $
;     postscript=postscript, atlascolor=intcolor2, atlaspsize=intpsize, atlassym=intsym

    djs_oplot, massaxis, lconst1, line=1, thick=postthick2
    djs_oplot, massaxis, lconst2, line=2, thick=postthick2
    djs_oplot, massaxis, lconst3, line=3, thick=postthick2

; ##########
; 0.2<z<0.4
; ##########

    bin1 = where((zobj gt 0.2) and (zobj le 0.4),nbin1)
    hi = where(oh[bin1] gt agesohcut,nhi,comp=lo,ncomp=nlo)

    x = mass[bin1]
    xerr = masserr[bin1]
    y = ml[bin1]
    yerr = mlerr[bin1]

    medyerr = median(yerr)
    medxerr = median(xerr)

; draw "lo" before "hi" because of the colors and the point density    
    
    ages_lineplot, x[lo], y[lo], xerr[lo], yerr[lo], plottype=1, postscript=postscript, $
      xtitle='', ytitle='', xrange=xrange, yrange=yrange, $
      position=pos[*,1], charsize=charsize_6, xtickname=replicate(' ',10), $
      ytickname=replicate(' ',10), /noerase, $
      agespsize=agesohlopsize, agessym=agesohlosym, agescolor=agesohlocolor
    if (nhi ne 0L) then ages_lineplot, x[hi], y[hi], xerr[hi], yerr[hi], $
      /overplot, agespsize=agesohhipsize, agessym=agesohhisym, agescolor=agesohhicolor
    legend, '0.2<z<0.4', /left, /top, box=0, charsize=charsize_2, charthick=postthick
;   oploterror, !x.crange[1]-yerrpos, !y.crange[0]+yerrpos*1.2, medxerr, $
;     medyerr, ps=3, /data, thick=postthick3, errthick=postthick

;   atlas1d_lineplot, massint, mlint, massinterr, mlinterr, /overplot, $
;     postscript=postscript, atlascolor=intcolor2, atlaspsize=intpsize, atlassym=intsym

    djs_oplot, massaxis, lconst1, line=1, thick=postthick2
    djs_oplot, massaxis, lconst2, line=2, thick=postthick2
    djs_oplot, massaxis, lconst3, line=3, thick=postthick2

; ##########
; 0.4<z<0.6
; ##########

    bin1 = where((zobj gt 0.4) and (zobj le 0.6),nbin1)
    hi = where(oh[bin1] gt agesohcut,nhi,comp=lo,ncomp=nlo)

    x = mass[bin1]
    xerr = masserr[bin1]
    y = ml[bin1]
    yerr = mlerr[bin1]

    medyerr = median(yerr)
    medxerr = median(xerr)

    ages_lineplot, x[hi], y[hi], xerr[hi], yerr[hi], plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, $
      position=pos[*,2], charsize=charsize_6, /noerase, $
      agespsize=agesohhipsize*1.4, agessym=agesohhisym, agescolor=agesohhicolor
    if (nlo ne 0L) then ages_lineplot, x[lo], y[lo], xerr[lo], yerr[lo], $
      /overplot, agespsize=agesohlopsize*1.2, agessym=agesohlosym, agescolor=agesohlocolor
    legend, '0.4<z<0.6', /left, /top, box=0, charsize=charsize_2, charthick=postthick
;   oploterror, !x.crange[1]-yerrpos, !y.crange[0]+yerrpos*1.2, medxerr, $
;     medyerr, ps=3, /data, thick=postthick3, errthick=postthick

;   atlas1d_lineplot, massint, mlint, massinterr, mlinterr, /overplot, $
;     postscript=postscript, atlascolor=intcolor2, atlaspsize=intpsize, atlassym=intsym

    djs_oplot, massaxis, lconst1, line=1, thick=postthick2
    djs_oplot, massaxis, lconst2, line=2, thick=postthick2
    djs_oplot, massaxis, lconst3, line=3, thick=postthick2

; ##########
; 0.6<z<0.8
; ##########

    bin1 = where((zobj gt 0.6) and (zobj le max(zobj)),nbin1)
    hi = where(oh[bin1] gt agesohcut,nhi,comp=lo,ncomp=nlo)

    x = mass[bin1]
    xerr = masserr[bin1]
    y = ml[bin1]
    yerr = mlerr[bin1]

    medyerr = median(yerr)
    medxerr = median(xerr)

    ages_lineplot, x[hi], y[hi], xerr[hi], yerr[hi], plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle='', xrange=xrange, yrange=yrange, $
      position=pos[*,3], charsize=charsize_6, /noerase, ytickname=replicate(' ',10), $
      agespsize=agesohhipsize*1.4, agessym=agesohhisym, agescolor=agesohhicolor
    if (nlo ne 0L) then ages_lineplot, x[lo], y[lo], xerr[lo], yerr[lo], $
      /overplot, agespsize=agesohlopsize*1.2, agessym=agesohlosym, agescolor=agesohlocolor
    legend, '0.6<z<0.8', /left, /top, box=0, charsize=charsize_2, charthick=postthick
;   oploterror, !x.crange[1]-yerrpos, !y.crange[0]+yerrpos*1.2, medxerr, $
;     medyerr, ps=3, /data, thick=postthick3, errthick=postthick

;   atlas1d_lineplot, massint, mlint, massinterr, mlinterr, /overplot, $
;     postscript=postscript, atlascolor=intcolor2, atlaspsize=intpsize, atlassym=intsym

    djs_oplot, massaxis, lconst1, line=1, thick=postthick2
    djs_oplot, massaxis, lconst2, line=2, thick=postthick2
    djs_oplot, massaxis, lconst3, line=3, thick=postthick2

    im_openclose, postscript=postscript, /close    

