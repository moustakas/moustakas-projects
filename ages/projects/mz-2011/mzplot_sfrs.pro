pro mzplot_sfrs, ps=ps
; jm11may17ucsd - SFR and SFR/M plots

    mzpath = mz_path()
    qapath = mzpath+'qaplots/'
    if keyword_set(ps) then begin
       pspath = qapath
       suffix = '.ps'
    endif else begin
       pspath = mz_path(/paper)
       suffix = '.eps'
    endelse

    sdsssfrs = read_mz_sample(/sfrs,/sdss)
    agessfrs = read_mz_sample(/sfrs)
    zbins = mz_zbins(nzbins)

; --------------------------------------------------
; QAplot: SDSS - SFR(MPA-DR7) vs SFR(MPA-DR4)
    psfile = qapath+'sfrmpa_dr7_vs_sfrmpa_dr4.ps'
    im_plotconfig, 6, pos, psfile=psfile, yspace=1.0, xmargin=[1.3,0.4], $
      width=6.8, height=[5.5,3.0]

    xrange = [-2.2,2.5]
    yrange = [-2.2,2.5]
    massrange = [8.0,12.0]
    residrange = 2.8*[-1,1]
    sfraxis = range(xrange[0],xrange[1],50)
    
    good = where((sdsssfrs.agn eq 0) and (sdsssfrs.sfrmpa_dr4_50 gt -900) and $
      (sdsssfrs.sfrmpa_flag eq 0),ngood)

    mass = sdsssfrs[good].mass_50
    xx = sdsssfrs[good].sfrmpa_50
    yy = sdsssfrs[good].sfrmpa_dr4_50
    xxerr = sdsssfrs[good].sfrmpa_err
    yyerr = sdsssfrs[good].sfrmpa_dr4_err
    resid = yy-xx
    splog, djs_median(resid), djsig(resid,sigrej=4.0)
    
    mzplot_scatterplot, xx, yy, /sdss, ysty=1, xsty=1, position=pos[*,0], $
      xrange=xrange, yrange=yrange, xtitle=mzplot_sfrtitle()+' [B04-DR7]', $
      ytitle=mzplot_sfrtitle()+' [B04-DR4]', levels=[0.5,0.9,0.975]
    oploterror, xrange[1]-median(xxerr)-0.5, yrange[0]+median(yyerr)+0.3, $
      median(xxerr), median(yyerr), errthick=6
    djs_oplot, sfraxis, sfraxis, line=0, thick=4

    mzplot_scatterplot, mass, resid, /sdss, /noerase, ysty=1, xsty=1, position=pos[*,1], $
      xrange=massrange, yrange=residrange, xtitle=mzplot_masstitle(), ytitle='Residuals (dex)', $
      levels=[0.5,0.9,0.975]
    djs_oplot, !x.crange, [0,0], line=0, thick=4
    
    im_plotconfig, /psclose, psfile=psfile, /gzip
;   im_plotconfig, /psclose, psfile=psfile, gzip=keyword_set(ps)

; --------------------------------------------------
; QAplot: SDSS - SFR(MPA) vs SFR(iSEDfit)
    psfile = qapath+'sdss_sfrmpa_vs_sfrisedfit.ps'
;   psfile = pspath+'sdss_sfrmpa_vs_sfrisedfit'+suffix
    im_plotconfig, 6, pos, psfile=psfile, yspace=1.0, xmargin=[1.3,0.4], $
      width=6.8, height=[5.5,3.0]

    xrange = [-2.2,2.2]
    yrange = [-2.2,2.2]
    massrange = [8.0,12.0]
    residrange = 2.8*[-1,1]
    sfraxis = range(xrange[0],xrange[1],50)
    
    good = where((sdsssfrs.agn eq 0) and (sdsssfrs.sfrmpa_flag eq 0),ngood)
;   good = where((sdsssfrs.agn eq 0) and (sdsssfrs.ohsample eq 1) and (sdsssfrs.sfrmpa_flag eq 0),ngood)

    mass = sdsssfrs[good].mass_50
;   xx = sdsssfrs[good].sfr_avg
;   xx = sdsssfrs[good].sfr100_avg
    xx = sdsssfrs[good].sfr100_50
    yy = sdsssfrs[good].sfrmpa_50
    xxerr = sdsssfrs[good].sfr100_err
    yyerr = sdsssfrs[good].sfrmpa_err
    resid = yy-xx
    med = im_medxbin(mass,resid,0.1,minpt=500,/verbose)
    splog, djs_median(resid), djsig(resid,sigrej=4.0)
    
    mzplot_scatterplot, xx, yy, /sdss, ysty=1, xsty=1, position=pos[*,0], $
      xrange=xrange, yrange=yrange, xtitle=mzplot_sfrtitle()+' [iSEDfit]', $
      ytitle=mzplot_sfrtitle()+' [B04-DR7]', levels=[0.5,0.9,0.975]
    oploterror, xrange[1]-median(xxerr)-0.5, yrange[0]+median(yyerr)+0.3, $
      median(xxerr), median(yyerr), errthick=6
    djs_oplot, sfraxis, sfraxis, line=0, thick=4

    mzplot_scatterplot, mass, resid, /sdss, /noerase, ysty=1, xsty=1, position=pos[*,1], $
      xrange=massrange, yrange=residrange, xtitle=mzplot_masstitle(), ytitle='Residuals (dex)', $
      levels=[0.5,0.9,0.975]
    oploterror, med.xbin, med.medy, med.sigy, psym=-6, color=djs_icolor('red'), errcolor=djs_icolor('red')
    djs_oplot, !x.crange, [0,0], line=0, thick=4
    
    im_plotconfig, /psclose, psfile=psfile, /gzip
;   im_plotconfig, /psclose, psfile=psfile, gzip=keyword_set(ps)

; --------------------------------------------------
; QAplot: AGES - SFR(24) vs SFR(iSEDfit)
    psfile = qapath+'ages_sfr24_vs_sfrisedfit.ps'
;   psfile = pspath+'ages_sfr24_vs_sfrisedfit'+suffix
    im_plotconfig, 4, pos, psfile=psfile, yspace=[1.0,1.0], xmargin=[1.3,0.4], $
      width=6.8, height=[3.0,2.0,2.0], charsize=1.5
;   im_plotconfig, 6, pos, psfile=psfile, yspace=1.0, xmargin=[1.3,0.4], $
;     width=6.8, height=[5.5,3.0]

    xrange = [-1.4,2.7]
    yrange = [-1.4,2.7]
    massrange = [9.0,12.0]
    zrange = [0.0,0.8]
    residrange = 2.3*[-1,1]
    sfraxis = range(xrange[0],xrange[1],50)

;   good = lindgen(n_elements(agessfrs))
    good = where(agessfrs.mips,ngood)
    mass = agessfrs[good].mass_50
    zobj = agessfrs[good].z
    weight = agessfrs[good].weight
    xx = agessfrs[good].sfr100_avg
;   xx = agessfrs[good].sfr100_50
    yy = agessfrs[good].sfr24
    xxerr = agessfrs[good].sfr100_err
    yyerr = agessfrs[good].sfr24_err
    resid = yy-xx
    splog, djs_median(resid), djsig(resid,sigrej=4.0)

    med_mass = im_medxbin(mass,resid,0.1,minpt=10,/verbose)
    med_zobj = im_medxbin(zobj,resid,0.03,minpt=10,/verbose)

    sf = where(agessfrs[good].agn eq 0)
    unk = where(agessfrs[good].agn eq -1)
    agn = where(agessfrs[good].agn eq 1)
    help, good, sf, unk, agn
    
;   mzplot_scatterplot, xx[sf], yy[sf], ysty=1, xsty=1, npix=30, position=pos[*,0], $
    mzplot_scatterplot, xx, yy, weight=weight, ysty=1, xsty=1, npix=30, position=pos[*,0], $
      xrange=xrange, yrange=yrange, xtitle=mzplot_sfrtitle()+' [iSEDfit]', $
      ytitle=mzplot_sfrtitle()+' [L(24)]', levels=[0.5,0.75,0.9]
;   djs_oplot, xx[agn], yy[agn], psym=6, sym=0.1, color='dark green'
;   djs_oplot, xx[unk], yy[unk], psym=6, sym=0.1, color='blue'
    oploterror, xrange[1]-median(xxerr)-0.3, yrange[0]+median(yyerr)+0.3, $
      median(xxerr), median(yyerr), errthick=6
    djs_oplot, sfraxis, sfraxis, line=0, thick=4

; residuals vs mass
    mzplot_scatterplot, mass, resid, weight=weight, /noerase, ysty=1, xsty=1, position=pos[*,1], $
;   mzplot_scatterplot, mass[sf], resid[sf], /noerase, ysty=1, xsty=1, position=pos[*,1], $
      xrange=massrange, yrange=residrange, xtitle=mzplot_masstitle(), $
      ytitle='Residuals (dex)', levels=[0.5,0.75,0.9], yminor=3
;   djs_oplot, mass[agn], resid[agn], psym=6, sym=0.1, color='dark green'
;   djs_oplot, mass[unk], resid[unk], psym=6, sym=0.1, color='blue'
    oploterror, med_mass.xbin, med_mass.medy, med_mass.sigy, psym=-6, $
      color=djs_icolor('red'), errcolor=djs_icolor('red')
    djs_oplot, !x.crange, [0,0], line=0, thick=4
    
; residuals vs redshift
    mzplot_scatterplot, zobj, resid, weight=weight, /noerase, ysty=1, xsty=1, position=pos[*,2], $
      xrange=zrange, yrange=residrange, xtitle='Redshift', $
      ytitle='Residuals (dex)', levels=[0.5,0.75,0.9], yminor=3
    oploterror, med_zobj.xbin, med_zobj.medy, med_zobj.sigy, psym=-6, $
      color=djs_icolor('red'), errcolor=djs_icolor('red')
    djs_oplot, !x.crange, [0,0], line=0, thick=4
    
    im_plotconfig, /psclose, psfile=psfile, /gzip
;   im_plotconfig, /psclose, psfile=psfile, gzip=keyword_set(ps)

; --------------------------------------------------
; paper plot: AGES - SFR(24)/M vs M in six bins of redshift
    psfile = pspath+'ages_sfrm_vs_mass'+suffix
    im_plotconfig, 7, pos, psfile=psfile, charsize=1.5, $
      height=2.7*[1,1,1]

    xrange = [9.1,11.9]
    yrange = [-2.4,+0.7]
    xtitle1 = mzplot_masstitle()
    ytitle1 = mzplot_sfrmtitle()
    maxis = range(9.0,12.5,50)

;   good = lindgen(n_elements(agessfrs))
    good = where(agessfrs.mips,ngood)
    zobj = agessfrs[good].z
    weight = agessfrs[good].weight
    mass = agessfrs[good].mass_50
    sfrm = (agessfrs[good].sfr24+9-mass)>(-4) ; [Gyr^-1]
    sfrm_ised = (agessfrs[good].sfr100_50+9-mass)>(-4) ; [Gyr^-1]

    ohgood = where(agessfrs.mips and agessfrs.ohsample)
    ohzobj = agessfrs[ohgood].z
    ohweight = agessfrs[ohgood].weight
    ohmass = agessfrs[ohgood].mass_50
    ohsfrm = (agessfrs[ohgood].sfr24+9-ohmass)>(-4) ; [Gyr^-1]
    ohsfrm_ised = (agessfrs[ohgood].sfr100_50+9-ohmass)>(-4) ; [Gyr^-1]
    
    for iz = 0, nzbins-1 do begin
       these = where((zobj ge zbins[iz].zlo) and $
         (zobj lt zbins[iz].zup),nthese)
       xx = mass[these]
       yy = sfrm[these]
       yy_ised = sfrm_ised[these]
       ww = weight[these]

       ohthese = where((ohzobj ge zbins[iz].zlo) and $
         (ohzobj lt zbins[iz].zup),nohthese)
       ohxx = ohmass[ohthese]
       ohyy = ohsfrm[ohthese]
       ohyy_ised = ohsfrm_ised[ohthese]
       ohww = ohweight[ohthese]

       if odd(iz) then begin
          ytitle = '' & ytickname = replicate(' ',10)
       endif else begin
          ytitle = ytitle1 & delvarx, ytickname
       endelse
       if (iz lt 4) then begin
          xtitle = '' & xtickname = replicate(' ',10)
       endif else begin
          xtitle = xtitle1
          delvarx, xtickname
       endelse
       
       mzplot_scatterplot, xx, yy, weight=ww, ysty=1, xsty=1, noerase=(iz gt 0), $
         position=pos[*,iz], xrange=xrange, yrange=yrange, xtitle=xtitle, $
         ytitle=ytitle, xtickname=xtickname, ytickname=ytickname, $
         levels=[0.5,0.75,0.9]
       mzplot_scatterplot, ohxx, ohyy, weight=ohww, ysty=1, xsty=1, /noerase, $
         position=pos[*,iz], xrange=xrange, yrange=yrange, xtitle='', $
         ytitle='', xtickname=replicate(' ',10), ytickname=replicate(' ',10), $
         levels=[0.5,0.75,0.9], /nogrey, ccolor=djs_icolor('blue'), /nooutlier, $
         npix=8
;      mzplot_scatterplot, xx, yy_ised, weight=ww, ysty=1, xsty=1, /noerase, $
;        position=pos[*,iz], xrange=xrange, yrange=yrange, xtitle='', $
;        ytitle='', xtickname=replicate(' ',10), ytickname=replicate(' ',10), $
;        levels=[0.5,0.75,0.9], /nogrey, ccolor=djs_icolor('red'), /nooutlier
       djs_oplot, maxis, poly(maxis,[-6.33,-0.35])+9, line=0, thick=5                ; Salim+07
       djs_oplot, maxis, poly(maxis,[-6.33,-0.35])+9+alog10(3.0), line=5, thick=5    ; factor of 3
       legend, 'z='+strtrim(string(zbins[iz].zlo,format='(F12.2)'),2)+$
         '-'+strtrim(string(zbins[iz].zup,format='(F12.2)'),2), $
         /right, /top, box=0, margin=0, charsize=1.4

       djs_oplot, [9.25,9.75,10.25,10.75,11.25], [-0.7,-0.8,-0.9,-1.1,-1.8], $
         psym=6, sym=3
       
    endfor
       
    im_plotconfig, /psclose, psfile=psfile, gzip=keyword_set(ps)

; --------------------------------------------------
; paper plot: AGES - SFR(24) vs M in six bins of redshift
    psfile = pspath+'ages_sfr_vs_mass'+suffix
    im_plotconfig, 7, pos, psfile=psfile, charsize=1.5, $
      height=2.7*[1,1,1]

    xrange = [9.3,11.9]
    yrange = [-1.3,3.2]
    xtitle1 = mzplot_masstitle()
    ytitle1 = mzplot_sfrtitle()
    maxis = range(xrange[0]+0.1,xrange[1]-0.1,50)

;   good = lindgen(n_elements(agessfrs))
    good = where(agessfrs.mips,ngood)
    zobj = agessfrs[good].z
    weight = agessfrs[good].weight
    mass = agessfrs[good].mass_50
    sfr = agessfrs[good].sfr24
    sfr_ised = agessfrs[good].sfr100_50

    ohgood = where(agessfrs.mips and agessfrs.ohsample)
    ohzobj = agessfrs[ohgood].z
    ohweight = agessfrs[ohgood].weight
    ohmass = agessfrs[ohgood].mass_50
    ohsfr = agessfrs[ohgood].sfr24
    ohsfr_ised = agessfrs[ohgood].sfr100_50
    
    for iz = 0, nzbins-1 do begin
       these = where((zobj ge zbins[iz].zlo) and $
         (zobj lt zbins[iz].zup),nthese)
       xx = mass[these]
       yy = sfr[these]
       yy_ised = sfr_ised[these]
       ww = weight[these]

       ohthese = where((ohzobj ge zbins[iz].zlo) and $
         (ohzobj lt zbins[iz].zup),nohthese)
       ohxx = ohmass[ohthese]
       ohyy = ohsfr[ohthese]
       ohyy_ised = ohsfr_ised[ohthese]
       ohww = ohweight[ohthese]

       if odd(iz) then begin
          ytitle = '' & ytickname = replicate(' ',10)
       endif else begin
          ytitle = ytitle1 & delvarx, ytickname
       endelse
       if (iz lt 4) then begin
          xtitle = '' & xtickname = replicate(' ',10)
       endif else begin
          xtitle = xtitle1
          delvarx, xtickname
       endelse
       
       mzplot_scatterplot, xx, yy, weight=ww, ysty=1, xsty=1, noerase=(iz gt 0), $
         position=pos[*,iz], xrange=xrange, yrange=yrange, xtitle=xtitle, $
         ytitle=ytitle, xtickname=xtickname, ytickname=ytickname, $
         levels=[0.5,0.75,0.9]
       mzplot_scatterplot, ohxx, ohyy, weight=ohww, ysty=1, xsty=1, /noerase, $
         position=pos[*,iz], xrange=xrange, yrange=yrange, xtitle='', $
         ytitle='', xtickname=replicate(' ',10), ytickname=replicate(' ',10), $
         levels=[0.5,0.75,0.9], /nogrey, ccolor=djs_icolor('blue'), /nooutlier, $
         npix=8
;      djs_oplot, maxis, poly(maxis,[-6.33,-0.35])+9, line=0, thick=5                ; Salim+07
;      djs_oplot, maxis, poly(maxis,[-6.33,-0.35])+9+alog10(3.0), line=5, thick=5    ; factor of 3
       legend, 'z='+strtrim(string(zbins[iz].zlo,format='(F12.2)'),2)+$
         '-'+strtrim(string(zbins[iz].zup,format='(F12.2)'),2), $
         /left, /top, box=0, margin=0, charsize=1.4
       djs_oplot, maxis, poly(maxis-10,[alog10(1.48),0.65]), line=0, thick=5
       djs_oplot, maxis, poly(maxis-10,[alog10(7.5),0.67]), line=5, thick=5

;      djs_oplot, [9.25,9.75,10.25,10.75,11.25], [-0.7,-0.8,-0.9,-1.1,-1.8], $
;        psym=6, sym=3
       
    endfor
       
    im_plotconfig, /psclose, psfile=psfile, gzip=keyword_set(ps)

; --------------------------------------------------
; paper plot: SDSS - SFR/M vs M
    psfile = pspath+'sdss_sfrm_vs_mass'+suffix
    im_plotconfig, 0, pos, psfile=psfile, charsize=2.0, $
      width=6.7, xmargin=[1.4,0.4], height=6.3

    xrange = [8.3,11.9]
    yrange = [-3.4,+0.8]
    xtitle = mzplot_masstitle()
    ytitle = mzplot_sfrmtitle()
    maxis = range(8.0,12.5,50)

    good = where((sdsssfrs.agn eq 0) and (sdsssfrs.sfrmpa_flag eq 0),ngood)
    mass = sdsssfrs[good].massmpa_50
;   mass = sdsssfrs[good].mass
;   sfrm = (sdsssfrs[good].sfr100_avg-mass+9)>(-4)
    sfrm = (sdsssfrs[good].sfrmpa_50-mass+9)>(-4)
    weight = sdsssfrs[good].weight
    
    mzplot_scatterplot, mass, sfrm, weight=weight, ysty=1, xsty=1, $
      /sdss, position=pos, xrange=xrange, yrange=yrange, xtitle=xtitle, $
      ytitle=ytitle;, levels=[0.5,0.75,0.9]
    djs_oplot, maxis, poly(maxis,[-6.33,-0.35])+9, line=0, thick=5 ; Salim+07 (Chabrier IMF)
       
    im_plotconfig, /psclose, psfile=psfile, gzip=keyword_set(ps)

return
end
