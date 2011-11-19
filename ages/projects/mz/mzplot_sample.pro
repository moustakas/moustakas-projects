pro ages_oplot, x, y, symsize=symsize, color=color, _extra=extra
    if (n_elements(symsize) eq 0) then symsize = 0.4
    if (n_elements(color) eq 0) then color = 'grey'
    djs_oplot, x, y, psym=symcat(16), symsize=symsize, $
      color=im_color(color), _extra=extra
return
end

pro mzplot_sample, ps=ps
; jm09mar25nyu - sample definition plots

    common com_quick, sdss_dfactor, taumodel
    
    mzpath = mz_path()
    qapath = mzpath+'qaplots/'
    if keyword_set(ps) then begin
       pspath = qapath
       suffix = '.ps'
    endif else begin
       pspath = mz_path(/paper)
       suffix = '.eps'
    endelse

    ocor = 1.0 + 1.0/2.984 ; intrinsic 5007/4959 ratio
    lsun = 3.826D33        ; [erg/s]
    dist = 3.085678D19     ; fiducial distance [10 pc in cm]

; ------------------------------------------------------------
; Figures 2 & 3 - AGES - H-beta, [OII], and [OIII] selection
    agesparent = read_mz_sample(/parent)
    agesmass = read_mz_sample(/mass)
    agesispec = read_ages(/ppxf)

    match, agesparent.ages_id, agesispec.ages_id, m1, m2
    srt = sort(m1) & m1 = m1[srt] & m2 = m2[srt]
    agesispec = agesispec[m2]

    ages_hbcut = mz_hbcut()
    ages_snrcut = 1.0
    ages_oiihbcut = -0.3
;   ages_ewcut = 0.0
;   ages_sigmacut = 400.0 

    zaxis = range(0.01,1.0,100)
    scale = 1D-40
    zrange = [0.0,0.8]
    ages_dfactor = 4.0*!dpi*dluminosity(agesispec.z,/cm)^2
    
    sel1 = where($
      (agesispec.h_beta[0] gt ages_hbcut) and $
      (agesispec.h_beta[0]/agesispec.h_beta[1] gt ages_snrcut) and $ ; nominal cut
;     (agesispec.h_beta_ew[0] gt ages_ewcut) and $
;     (agesispec.h_beta_sigma[0] lt ages_sigmacut) and $ ; nominal cut
      (agesispec.oii_3727[1] ne -2.0))
    rej1 = where($
      (agesispec.h_beta[0] le ages_hbcut) and $
      (agesispec.h_beta[0]/agesispec.h_beta[1] gt ages_snrcut) and $ ; nominal cut
;     (agesispec.h_beta_ew[0] gt ages_ewcut) and $
;     (agesispec.h_beta_sigma[0] lt ages_sigmacut) and $
      (agesispec.oii_3727[1] ne -2.0))

    sel2 = where($
      (agesispec.h_beta[0] gt ages_hbcut) and $
;     (agesispec.h_beta_ew[0] gt ages_ewcut) and $
;     (agesispec.h_beta_sigma[0] lt ages_sigmacut) and $
      (agesispec.oii_3727[1] ne -2.0) and $
      (agesispec.oii_3727[0] gt 10^ages_oiihbcut*agesispec.h_beta[0]))
    rej2 = where($
      (agesispec.h_beta[0] gt ages_hbcut) and $
;     (agesispec.h_beta_ew[0] gt ages_ewcut) and $
;     (agesispec.h_beta_sigma[0] lt ages_sigmacut) and $
      (agesispec.oii_3727[1] ne -2.0) and $
      (agesispec.oii_3727[0] le 10^ages_oiihbcut*agesispec.h_beta[0]) and $
      (agesispec.oii_3727[1] gt 0.0))
    help, sel1, sel2
    
; Figure 2 - redshift vs H-beta
    psfile = pspath+'z_vs_hb'+suffix
    im_plotconfig, 6, pos, psfile=psfile, xmargin=[1.7,0.4], $
      width=6.0, height=[3.0,3.0]

; luminosity    
    djs_plot, [0], [0], /nodata, position=pos[*,0], xsty=1, ysty=1, $
      xrange=zrange, yrange=scale*[1.05D38,1.5D42], xtitle='', xtickname=replicate(' ',10), $
      ytitle='L(H\beta) (10^{40} erg s^{-1})', /ylog
    djs_oplot, agesispec[sel1].z, scale*ages_dfactor[sel1]*agesispec[sel1].h_beta[0], $
      psym=symcat(9,thick=1), symsize=0.3, color=im_color('grey60')
    djs_oplot, agesispec[rej1].z, scale*ages_dfactor[rej1]*agesispec[rej1].h_beta[0], $
      psym=symcat(16), symsize=0.6, color=im_color('midnight blue')
    djs_oplot, zaxis, scale*ages_hbcut*4.0*!dpi*dluminosity(zaxis,/cm)^2, line=0, thick=6
    im_legend, 'F(H\beta)>3\times10^{-17} erg s^{-1} cm^{-2}', $
      /right, /bottom, box=0, charsize=1.4, margin=0, line=0, thick=6

; EW
    djs_plot, [0], [0], /nodata, /noerase, position=pos[*,1], xsty=1, ysty=1, $
      xrange=zrange, yrange=[0.1,150], xtitle='Redshift', $
      ytitle='EW(H\beta) (\AA)', /ylog
    djs_oplot, agesispec[sel1].z, agesispec[sel1].h_beta_ew[0], $
      psym=symcat(9,thick=1), symsize=0.3, color=im_color('grey60')
    djs_oplot, agesispec[rej1].z, agesispec[rej1].h_beta_ew[0], $
      psym=symcat(16), symsize=0.6, color=im_color('midnight blue')
;   djs_oplot, zaxis, ages_hbcut*4.0*!dpi*dluminosity(zaxis,/cm)^2/ages_l4861_limit, line=0, thick=6
    im_plotconfig, /psclose, psfile=psfile

; Figure 3 - [OII]/Hb and [OIII]/Hb vs redshift
    psfile = pspath+'z_vs_oiihb_oiiihb'+suffix
    im_plotconfig, 6, pos, psfile=psfile, xmargin=[1.7,0.4], $
      width=6.0, height=[3.0,3.0], charsize=1.6

; [OII]/Hb    
    djs_plot, [0], [0], /nodata, position=pos[*,0], xsty=1, ysty=1, $
      xrange=zrange, yrange=[-1.5,1.0], xtitle='', xtickname=replicate(' ',10), $
      ytitle='log ([O II]/H\beta)'

    lim = where(agesispec[sel1].oii_3727[1] eq -1.0,comp=det)
    djs_oplot, agesispec[sel1[det]].z, alog10(agesispec[sel1[det]].oii_3727[0]/agesispec[sel1[det]].h_beta[0]), $
      psym=symcat(9,thick=1), symsize=0.3, color=im_color('grey60')

    plotsym, 1.0, 1.5, thick=4
    djs_oplot, agesispec[sel1[lim]].z, alog10(agesispec[sel1[lim]].oii_3727_limit/agesispec[sel1[lim]].h_beta[0]), $
      psym=8, color=im_color('orange red',101)
    djs_oplot, !x.crange, ages_oiihbcut*[1,1], line=5, thick=6
    
; [OIII]/Hb
    djs_plot, [0], [0], /nodata, /noerase, position=pos[*,1], xsty=1, ysty=1, $
      xrange=zrange, yrange=[-1.8,1.5], xtitle='Redshift', $
      ytitle='log ([O III] \lambda5007/H\beta)'

    lim = where(agesispec[sel2].oiii_5007[1] eq -1.0,comp=det)
    djs_oplot, agesispec[sel2[det]].z, alog10(agesispec[sel2[det]].oiii_5007[0]/$
      agesispec[sel2[det]].h_beta[0]), psym=symcat(9,thick=1), symsize=0.3, color=im_color('grey60')

    plotsym, 1.0, 0.8, thick=4
    djs_oplot, agesispec[sel2[lim]].z, alog10(agesispec[sel2[lim]].oiii_5007_limit/$
      agesispec[sel2[lim]].h_beta[0]), psym=8, color=im_color('orange red',101)
    
    im_plotconfig, /psclose, psfile=psfile

stop    
    
; ------------------------------------------------------------
; Figure 6 - distribution of global properties - AGES & SDSS

;   sdssparent = read_mz_sample(/parent,/sdss)
    sdsskcorr = read_mz_sample(/mzhii_ancillary,/sdss)
    sdssmass = read_mz_sample(/mzhii_mass,/sdss)

    agesparent = read_mz_sample(/parent)
    agesparentmass = read_mz_sample(/mass)
    ageskcorr = read_mz_sample(/mzhii_ancillary)
    agesmass = read_mz_sample(/mzhii_mass)

    psfile = pspath+'ages_histograms'+suffix
    im_plotconfig, 5, pos, psfile=psfile, height=[3.0,3.0], $
      width=[3.0,3.0], xmargin=[1.2,1.2], xspace=0.1, yspace=1.0, $
      charsize=1.7

    zbinsize = 0.02 ; 0.01
    absmagbinsize = 0.2 ; 0.25
    colorbinsize = 0.03 ; 0.04
    massbinsize = 0.1 ; 0.15
    yfactor = 1.1

    splog, weighted_quantile(ageskcorr.z,ageskcorr.final_weight,quant=0.5), $
      im_weighted_mean(ageskcorr.z,weight=ageskcorr.final_weight)
    splog, weighted_quantile(sdsskcorr.z,sdsskcorr.final_weight,quant=0.5), $
      im_weighted_mean(sdsskcorr.z,weight=sdsskcorr.final_weight)

    histthick1 = 8
    fspacing = 0.04
    snorm = 15.0
;   snorm = n_elements(sdsskcorr.z)/float(n_elements(ageskcorr.z))                                                                

    scolor = 'dark red'
    sline = 0

    fullcolor = 'gray40'
    agescolor = 'powder blue'
    
; ####################
; Redshift
; ####################

    im_plothist, agesparent.z, xbin, ybin, bin=zbinsize, $
      weight=agesparent.final_weight, /noplot

    xtitle = 'Redshift'
    ytitle = 'Number'
    xrange = [0.0,0.8]
    yrange = [0,max(ybin)*yfactor]
    
    djs_plot, [0], [0], /nodata, xsty=5, ysty=5, xrange=xrange, yrange=yrange, $
      position=pos[*,0], xtickname=replicate(' ',10), ytickname=replicate(' ',10)
    im_plothist, agesparent.z, bin=zbinsize, weight=agesparent.final_weight, $
      /overplot, fcolor=im_color(fullcolor), color=im_color(fullcolor), thick=histthick1;, /fill, $
;     /fline, forientation=45, fspacing=fspacing
;     /overplot, /fill, fcolor=im_color(fullcolor), /fline, $
;     forientation=45, fspacing=fspacing, color=im_color(fullcolor)
    im_plothist, ageskcorr.z, bin=zbinsize, weight=ageskcorr.final_weight, $
      /overplot, thick=histthick1, /fill, color=im_color(agescolor), fcolor=im_color(agescolor)
    im_plothist, ageskcorr.z, bin=zbinsize, weight=ageskcorr.final_weight, /overplot, thick=5
    im_plothist, sdsskcorr.z, weight=sdsskcorr.final_weight, bin=zbinsize, /overplot, thick=histthick1, $
      line=sline, color=im_color(scolor), norm=snorm;, $
;     /fline, forientation=45, fspacing=0.1, /fill, fcolor=im_color(scolor)
    djs_plot, [0], [0], /nodata, /noerase, xtitle=xtitle, ytitle=ytitle, $
      xsty=1, ysty=1, xrange=xrange, yrange=yrange, position=pos[*,0], $
      ytickinterval=500
;   im_legend, '(a)', /left, /top, box=0, margin=0

; ####################
; MB
; ####################
    
    im_plothist, agesparent.k_ubvrijhk_absmag_00[1], bin=absmagbinsize, xbin, ybin, $
      weight=agesparent.final_weight, /noplot

    xtitle = mzplot_mbtitle()
    ytitle = 'Number'
    xrange = reverse([-16.1,-23.9])
    yrange = [0,max(ybin)*yfactor]
    
    djs_plot, [0], [0], /nodata, /noerase, xsty=5, ysty=5, xrange=xrange, yrange=yrange, $
      position=pos[*,1], xtickname=replicate(' ',10), ytickname=replicate(' ',10)
    im_plothist, agesparent.k_ubvrijhk_absmag_00[1], bin=absmagbinsize, weight=agesparent.final_weight, $
      /overplot, fcolor=im_color(fullcolor), color=im_color(fullcolor), thick=histthick1;, /fill, $
;     /fline, forientation=45, fspacing=fspacing
;     /overplot, /fill, fcolor=im_color(fullcolor), /fline, $
;     forientation=45, fspacing=fspacing, color=im_color(fullcolor)
    im_plothist, ageskcorr.k_ubvrijhk_absmag_00[1], bin=absmagbinsize, $
      weight=ageskcorr.final_weight, /overplot, thick=histthick1, /fill, color=im_color(agescolor), fcolor=im_color(agescolor)
    im_plothist, ageskcorr.k_ubvrijhk_absmag_00[1], bin=absmagbinsize, weight=ageskcorr.final_weight, /overplot, thick=5
    im_plothist, sdsskcorr.k_ubvrijhk_absmag_00[1], weight=sdsskcorr.final_weight, bin=absmagbinsize, $
      /overplot, thick=histthick1, line=sline, color=im_color(scolor,101), norm=snorm
    djs_plot, [0], [0], /nodata, /noerase, xtitle=xtitle, ytitle='', xsty=1, $
      ysty=1, xrange=xrange, yrange=yrange, position=pos[*,1], $
      ytickname=replicate(' ',10)
    axis, /yaxis, yrange=yrange, ystyle=1, ytitle=textoidl(ytitle), $
      ytickinterval=500
;   im_legend, '(b)', /left, /top, box=0, margin=0

; ####################
; (U-B) color
; ####################

    ub_parent = agesparent.k_ubvrijhk_absmag_00[0]-agesparent.k_ubvrijhk_absmag_00[1]
    ub = ageskcorr.k_ubvrijhk_absmag_00[0]-ageskcorr.k_ubvrijhk_absmag_00[1]
    sdss_ub = sdsskcorr.k_ubvrijhk_absmag_00[0]-sdsskcorr.k_ubvrijhk_absmag_00[1]

    im_plothist, ub_parent, xbin, ybin, bin=colorbinsize, $
      weight=agesparent.final_weight, /noplot

    xtitle = 'U - B'
    ytitle = 'Number'
    xrange = [0.45,1.55]
    yrange = [0,max(ybin)*yfactor]
    
    djs_plot, [0], [0], /nodata, /noerase, xsty=5, ysty=5, xrange=xrange, yrange=yrange, $
      position=pos[*,2], xtickname=replicate(' ',10), ytickname=replicate(' ',10), $
      xtickinterval=0.5
    im_plothist, ub_parent, bin=colorbinsize, weight=agesparent.final_weight, $
      /overplot, fcolor=im_color(fullcolor), color=im_color(fullcolor), thick=histthick1;, /fill, $
;     /fline, forientation=45, fspacing=fspacing
    im_plothist, ub, bin=colorbinsize, weight=ageskcorr.final_weight, /overplot, $
      thick=histthick1, /fill, color=im_color(agescolor), fcolor=im_color(agescolor)
    im_plothist, ub, bin=colorbinsize, weight=ageskcorr.final_weight, /overplot, thick=5
    im_plothist, sdss_ub, weight=sdsskcorr.final_weight, bin=colorbinsize, /overplot, thick=histthick1, $
      line=sline, color=im_color(scolor,101), norm=snorm
    djs_plot, [0], [0], /nodata, /noerase, xtitle=xtitle, xrange=xrange, ytitle=ytitle, $
      xsty=1, ysty=1, yrange=yrange, position=pos[*,2], ytickinterval=500, xtickinterval=0.2
;   im_legend, '(c)', /left, /top, box=0, margin=0

; ####################
; Stellar Mass
; ####################

    im_plothist, agesparentmass.mass_50, xbin, ybin, bin=massbinsize, $
      weight=agesparent.final_weight, /noplot

    xtitle = mzplot_masstitle()
    ytitle = 'Number'
    xrange = [8.1,11.9]
    yrange = [0,max(ybin)*yfactor]
    
    djs_plot, [0], [0], /nodata, /noerase, xsty=5, ysty=5, xrange=xrange, yrange=yrange, $
      position=pos[*,3], xtickname=replicate(' ',10), ytickname=replicate(' ',10)
    im_plothist, agesparentmass.mass_50, bin=massbinsize, weight=agesparent.final_weight, $
      /overplot, fcolor=im_color(fullcolor), color=im_color(fullcolor), thick=histthick1;, /fill, $
;     /fline, forientation=45, fspacing=fspacing
;     /overplot, /fill, fcolor=im_color(fullcolor), /fline, $
;     forientation=45, fspacing=fspacing, color=im_color(fullcolor)
    im_plothist, agesmass.mass_50, bin=massbinsize, /overplot, $
      weight=ageskcorr.final_weight, thick=histthick1, /fill, color=im_color(agescolor), fcolor=im_color(agescolor)
    im_plothist, agesmass.mass_50, bin=massbinsize, weight=ageskcorr.final_weight, /overplot, thick=5
    im_plothist, sdssmass.mass_50, weight=sdsskcorr.final_weight, bin=massbinsize, /overplot, $
      thick=histthick1, line=sline, color=im_color(scolor,101), norm=snorm
    djs_plot, [0], [0], /nodata, /noerase, xtitle=xtitle, ytitle='', $
      xsty=1, ysty=1, xrange=xrange, yrange=yrange, $
      position=pos[*,3], yminor=3, ytickname=replicate(' ',10)
    axis, /yaxis, yrange=yrange, ystyle=1, ytitle=textoidl(ytitle), $
      ytickinterval=500
;   im_legend, '(d)', /left, /top, box=0, margin=0

    im_plotconfig, /psclose, psfile=psfile

stop    
    
; ------------------------------------------------------------
; Figure 7 - redshift vs absolute magnitude and stellar mass
    ageskcorr = read_mz_sample(/mzhii_ancillary)
    agesmass = read_mz_sample(/mzhii_mass)

;; compute the tau model
;    if (n_elements(taumodel) eq 0) then begin
;       tau = 3.0
;       ssp = im_read_bc03(bc03_extras=ext,/silent)
;       tauflux = im_convolve_sfh(ssp,tau=tau,mstar=ext.m_,cspmstar=taumstar)
;       tauflux = lsun*tauflux/(4.0*!dpi*dist^2.0)/rebin(reform(taumstar,1,220),6900,220) ; [erg/s/cm2/A/Msun]
;       taumodel = {tau: tau, age: ssp.age, wave: ssp.wave, flux: tauflux, mstar: taumstar}
;    endif
;
;; interpolate    
;    nz = 100 
;    zz = range(0.02,1.0,nz)
;    tauzform = 2.0
;
;    tauindx = findex(taumodel.age/1D9,getage(zz)-getage(tauzform))
;    tauflux = interpolate(taumodel.flux,tauindx)
;    taumstar = interpolate(taumodel.mstar,tauindx)
;  
;    mbfilter = 'bessell_B.par'
;    ifaint = mz_ifaint(select_filter=ifilter)
;    maggies = rebin(reform(10^(-0.4*ifaint),1,1),1,nz)
;
;    model = replicate({z: 0.0, tau_MB: 0.0, tau_mstar: 0.0},nz)
;    model.z = zz
;    
;    kk = im_simple_kcorrect(zz,maggies,maggies*0.0+1.0,$
;      ifilter,mbfilter,taumodel.wave,tauflux,absmag=abs,$
;      scale=tauscale)
;    model.tau_mstar = alog10(taumstar*(dluminosity(zz,/cm)/3.085678D19)^2*tauscale)
;    model.tau_MB = reform(abs)

    zaxis = range(0.03,0.75,50)
    zbins = mz_zbins(nzbins)
    limits = mrdfits(mzpath+'mz_limits.fits.gz',1)
    
; now make the plot    
    psfile = pspath+'z_vs_mb_mass'+suffix
    im_plotconfig, 6, pos, psfile=psfile, height=3.0*[1,1], xmargin=[1.3,0.2], $
      ymargin=[0.8,1.0], charsize=1.9, width=5.5

    xrange = [0.04,0.76]
    yrange1 = [-15.3,-23.4]
    yrange2 = [7.7,11.9]

; redshift vs absolute magnitude
    plot, [0], [0], /nodata, ysty=1, xsty=1, xrange=xrange, yrange=yrange1, $
      xtitle='', ytitle=mzplot_mbtitle(), xtickname=replicate(' ',10), $
      position=pos[*,0], yminor=4, xtickinterval=0.1
    ages_oplot, ageskcorr.z, ageskcorr.k_ubvrijhk_absmag_00[1], color='dodger blue'
;   djs_oplot, model.z, model.tau_MB, line=0, $
;     color=im_color('firebrick',101), thick=8
;   djs_oplot, zaxis, poly(zaxis,limits.mblimit_50_coeff), line=0, thick=8
    for iz = 0, nzbins-1 do begin
       djs_oplot, [zbins[iz].zlo,zbins[iz].zup], limits.mblimit_50[iz]*[1,1], $
         line=0, thick=6
       djs_oplot, zbins[iz].zlo*[1,1], [limits.mblimit_50[iz],!y.crange[1]], $
         line=0, thick=6
       djs_oplot, zbins[iz].zup*[1,1], [limits.mblimit_50[iz],!y.crange[1]], $
         line=0, thick=6
    endfor
    
; redshift vs stellar mass
    plot, [0], [0], /nodata, /noerase, ysty=1, xsty=9, $
      xrange=xrange, yrange=yrange2, xtitle='Redshift', $
      ytitle=mzplot_masstitle(), position=pos[*,1], yminor=4, xtickinterval=0.1
    ages_oplot, ageskcorr.z, agesmass.mass_50, color='dodger blue'
    for iz = 0, nzbins-1 do begin
       djs_oplot, [zbins[iz].zlo,zbins[iz].zup], limits.masslimit_50[iz]*[1,1], $
         line=0, thick=6
       djs_oplot, zbins[iz].zlo*[1,1], [limits.masslimit_50[iz],!y.crange[1]], $
         line=0, thick=6
       djs_oplot, zbins[iz].zup*[1,1], [limits.masslimit_50[iz],!y.crange[1]], $
         line=0, thick=6
    endfor

;   djs_oplot, zaxis, poly(zaxis,limits.masslimit_50_coeff), line=0, thick=8
;   djs_oplot, model.z, model.tau_mstar, line=0, $
;     color=im_color('firebrick',101), thick=8
;   im_legend, ['\tau = 3 Gyr'], /right, /bottom, box=0, $
;     pspacing=1.4, color='firebrick', line=0, charsize=1.5, thick=8

    im_plotconfig, /psclose, psfile=psfile

; ------------------------------------------------------------
; Figure 1 - EW(Hb), EW([OII]), EW([O III]) and EW(R23) measured from
; our fluxed and unfluxed spectra
    agesparent = read_mz_sample(/parent)
    allispec = read_ages(/ppxf)
    unfluxed = read_ages(/ppxf,/unfluxed)

    match, agesparent.ages_id, allispec.ages_id, m1, m2
    srt = sort(m1) & m1 = m1[srt] & m2 = m2[srt]
    allispec = allispec[m2]
    unfluxed = unfluxed[m2]

    logewaxis = findgen((alog10(1E3)-alog10(1E-5))/0.05)*0.05+alog10(1E-5)
    ewaxis = 10.0^logewaxis
    
    psfile = pspath+'ages_ew_tests'+suffix
    im_plotconfig, 5, pos, psfile=psfile, height=2.95*[1,1], width=2.95*[1,1], $
      xmargin=[1.25,0.25], ymargin=[0.4,1.1], xspace=1.1, yspace=1.0, $
      charsize=1.5

    indx = where((allispec.h_beta_ew[0] gt 0.0) and (unfluxed.h_beta_ew[0] gt 0.0) and $
      (allispec.oii_3727_ew[0] gt 0.0) and (unfluxed.oii_3727_ew[0] gt 0.0) and $
      (allispec.oiii_5007_ew[0] gt 0.0) and (unfluxed.oiii_5007_ew[0] gt 0.0))
    
    z = allispec[indx].z
;   ewha = allispec[indx].h_alpha_ew[0]
    ewhb = allispec[indx].h_beta_ew[0]
    ewoii = allispec[indx].oii_3727_ew[0]
    ewoiii = allispec[indx].oiii_5007_ew[0]
    ewr23 = (ewoii+ocor*ewoiii)/ewhb

    ewoiicor = 1.2
    ewhbcor = 2.8
    dewhbcor = 1.5
;   ewhacor = 0.8
    ewoiiicor = 0.3
    
;   uewha = unfluxed[indx].h_alpha_ew[0]
    uewhb = unfluxed[indx].h_beta_ew[0]
    uewoii = unfluxed[indx].oii_3727_ew[0]
    uewoiii = unfluxed[indx].oiii_5007_ew[0]
    uewr23 = ((uewoii+ewoiicor)+ocor*(uewoiii+ewoiiicor))/(uewhb+ewhbcor)

; stats
;   oiistats = im_stats(100*(ewoii-uewoii)/(0.5*(ewoii+uewoii)),sigrej=3.0)
;   oiiistats = im_stats(100*(ewoiii-uewoiii)/(0.5*(ewoiii+uewoiii)),sigrej=3.0)
;   hbstats = im_stats(100*(ewhb-uewhb)/(0.5*(ewhb+uewhb)),sigrej=3.0)
;   hastats = im_stats(100*(ewha-uewha)/(0.5*(ewha+uewha)),sigrej=3.0)
;   r23stats = im_stats(100*(ewr23-uewr23)/(0.5*(ewr23+uewr23)),sigrej=3.0)
    oiistats = im_stats(alog10(ewoii/uewoii),sigrej=3.0)
    oiiistats = im_stats(alog10(ewoiii/uewoiii),sigrej=3.0)
    hbstats = im_stats(alog10(ewhb/uewhb),sigrej=3.0)
;   hastats = im_stats(alog10(ewha/uewha),sigrej=3.0)
    r23stats = im_stats(alog10(ewr23/uewr23),sigrej=3.0)

    splog, '[OII]:  '+$
      strtrim(string(oiistats.median_rej,format='(F12.4)'),2)+' ('+$
      strtrim(string(oiistats.mean_rej,format='(F12.4)'),2)+'+/-'+$
      strtrim(string(oiistats.sigma_rej,format='(F12.4)'),2)+')'
    splog, '[OIII]: '+$
      strtrim(string(oiiistats.median_rej,format='(F12.4)'),2)+' ('+$
      strtrim(string(oiiistats.mean_rej,format='(F12.4)'),2)+'+/-'+$
      strtrim(string(oiiistats.sigma_rej,format='(F12.4)'),2)+')'
    splog, 'Hb    : '+$
      strtrim(string(hbstats.median_rej,format='(F12.4)'),2)+' ('+$
      strtrim(string(hbstats.mean_rej,format='(F12.4)'),2)+'+/-'+$
      strtrim(string(hbstats.sigma_rej,format='(F12.4)'),2)+')'
;   splog, 'Ha    : '+$
;     strtrim(string(hastats.median_rej,format='(F12.4)'),2)+' ('+$
;     strtrim(string(hastats.mean_rej,format='(F12.4)'),2)+'+/-'+$
;     strtrim(string(hastats.sigma_rej,format='(F12.4)'),2)+')'
    splog, 'R23   : '+$
      strtrim(string(r23stats.median_rej,format='(F12.4)'),2)+' ('+$
      strtrim(string(r23stats.mean_rej,format='(F12.4)'),2)+'+/-'+$
      strtrim(string(r23stats.sigma_rej,format='(F12.4)'),2)+')'
    
; make the plot    
;   levels = errorf(0.5*(dindgen(2)+2)) ; 84.3%,96.7%
    levels = [0.5,0.95]
    
; EW(Hb)    
    xtitle = textoidl('EW(H\beta) (\AA) [Unfluxed]')
    ytitle = textoidl('EW(H\beta) (\AA) [Fluxed]')
    xrange = [0.15,99.0] & yrange = xrange
    xlogrange = alog10(xrange) & ylogrange = xlogrange

    mzplot_scatterplot, /ages, alog10(uewhb), alog10(ewhb), position=pos[*,0], $
      xtitle=xtitle, ytitle=ytitle, xstyle=1, ystyle=1, xrange=xlogrange, yrange=ylogrange, $
      levels=levels, ynpix=40, xnpix=40
    djs_oplot, logewaxis, logewaxis, line=0, thick=6.0
    djs_oplot, logewaxis, alog10(ewaxis+ewhbcor), line=0, color='dark red', thick=5
    djs_oplot, logewaxis, alog10(ewaxis+ewhbcor+dewhbcor), line=5, color='dark red', thick=5
    djs_oplot, logewaxis, alog10(ewaxis+ewhbcor-dewhbcor), line=5, color='dark red', thick=5
;   im_legend, '(a)', /left, /top, box=0, charsize=1.4, margin=0
    legend, textoidl('2.8\pm1.5 \AA'), /right, /bottom, box=0, charsize=1.4, $
      margin=0, line=0, thick=8.0, color=im_color('dark red'), pspacing=1.5
;   im_legend, '(a) H\beta', /left, /top, box=0, charsize=1.4, margin=0

;   mzages_hogg_scatterplot, z, uewhb-ewhb, outpsym=symcat(agessym2), $
;     outsymsize=agespsize2, outcolor=im_color(agescolor2,1E8), $
;     xtitle='Redshift', ytitle=textoidl('\Delta'+'EW (\AA)'), $
;     charsize=charsize_1, xsty=1, ysty=1, xrange=xlogrange, yrange=ylogrange, $
;     /internal, ynpix=30, xnpix=30, /outliers, /noerase, $
;     position=[pos[2,0]-0.13,pos[1,0]+0.06,pos[2,0]-0.01,pos[1,0]+0.12]
    
; EW([OIII])    
    xtitle = textoidl('EW([O III]) (\AA) [Unfluxed]')
    ytitle = textoidl('EW([O III]) (\AA) [Fluxed]')
    xrange = [0.3,220.0] & yrange = xrange
    xlogrange = alog10(xrange) & ylogrange = xlogrange

    mzplot_scatterplot, /ages, alog10(uewoiii), alog10(ewoiii), /noerase, position=pos[*,1], $
      xtitle=xtitle, ytitle=ytitle, xstyle=1, ystyle=1, xrange=xlogrange, yrange=ylogrange, $
      levels=levels, ynpix=40, xnpix=40
    djs_oplot, logewaxis, logewaxis, line=0, thick=6.0
;   djs_oplot, logewaxis, alog10(ewaxis+ewoiiicor), line=5, color='dark red', thick=8.0
;   im_legend, '(b)', /left, /top, box=0, charsize=1.4, margin=0
;   legend, textoidl('0.3 \AA'), /right, /bottom, box=0, charsize=1.4, $
;     margin=0, line=5, thick=8.0, color=im_color('dark red'), pspacing=1.5
;   im_legend, '(b) [O III] \lambda5007', /left, /top, box=0, charsize=1.4, margin=0

; EW([OII])    
    xtitle = textoidl('EW([O II]) (\AA) [Unfluxed]')
    ytitle = textoidl('EW([O II]) (\AA) [Fluxed]')
    xrange = [3.0,140.0] & yrange = xrange
    xlogrange = alog10(xrange) & ylogrange = xlogrange

    mzplot_scatterplot, /ages, alog10(uewoii), alog10(ewoii), /noerase, position=pos[*,2], $
      xtitle=xtitle, ytitle=ytitle, xstyle=1, ystyle=1, xrange=xlogrange, yrange=ylogrange, $
      levels=levels, ynpix=40, xnpix=40
    djs_oplot, logewaxis, logewaxis, line=0, thick=6.0
;   djs_oplot, logewaxis, alog10(ewaxis+ewoiicor), line=5, color='dark red', thick=linethick1
;   im_legend, '(c)', /left, /top, box=0, charsize=1.4, margin=0
;   im_legend, '(c) [O II]', /left, /top, box=0, charsize=1.4, margin=0

; EW(R23)
    xtitle = textoidl('log EW(R_{23}) [Unfluxed]')
    ytitle = textoidl('log EW(R_{23}) [Fluxed]')

    xrange = [0.6,23.0] & yrange = xrange
    xlogrange = alog10(xrange) & ylogrange = xlogrange

    mzplot_scatterplot, /ages, alog10(uewr23), alog10(ewr23), /noerase, position=pos[*,3], $
      xtitle=xtitle, ytitle=ytitle, xstyle=1, ystyle=1, xrange=xlogrange, yrange=ylogrange, $
      levels=levels, ynpix=40, xnpix=40
    djs_oplot, logewaxis, logewaxis, line=0, thick=6.0
;   im_legend, '(d)', /left, /top, box=0, charsize=1.4, margin=0
;   im_legend, '(d) R_{23}', /left, /top, box=0, charsize=1.4, margin=0
    
    im_plotconfig, /psclose, psfile=psfile

; ###########################################################################
; additional figures and QAplots

; ------------------------------------------------------------
; SDSS - H-beta, [OII], and [OIII] selection

    vagc = mz_get_vagc(sample=sample,letter=letter,poststr=poststr)
    sdssispec = read_vagc_garching(sample=sample,$
      letter=letter,poststr=poststr,/ispec)
    sdssparent = read_mz_sample(/parent,/sdss)
    keep = where(sdssispec.z gt 0.0)
    sdssispec = sdssispec[keep]
    sdssparent = sdssparent[keep]

; to figure out the limiting EW coefficients do:
    sdss_hbcut = mz_hbcut(/sdss)
    sdss_snrcut = 1.0
    sdss_oiihbcut = -0.3
;   sdss_ewcut = 0.0
;   sdss_sigmacut = 400.0

    zaxis = range(0.01,0.3,50)
    scale = 1D-40
    zrange = [0.0,0.25]
    if (n_elements(sdss_dfactor) eq 0L) then sdss_dfactor = 4.0*!dpi*dluminosity(sdssispec.z,/cm)^2

    sel1 = where($
      (sdssispec.h_beta[0] gt sdss_hbcut) and $
      (sdssispec.h_beta[0]/sdssispec.h_beta[1] gt sdss_snrcut) and $ ; nominal cut
;     (sdssispec.h_beta_ew[0] gt sdss_ewcut) and $
;     (sdssispec.h_beta_sigma[0] lt sdss_sigmacut) and $ ; nominal cut
      (sdssispec.oii_3727[1] ne -2.0) and $
      (sdssispec.oiii_5007[1] ne -2.0) and $ ; nominal
      (sdssispec.nii_6584[1] ne -2.0) and $ ; nominal
      (sdssispec.h_alpha[1] ne -2.0)) ; nominal
    rej1 = where($
      (sdssispec.h_beta[0] le sdss_hbcut) and $
      (sdssispec.h_beta[0]/sdssispec.h_beta[1] gt sdss_snrcut) and $ ; nominal cut
;     (sdssispec.h_beta_ew[0] gt sdss_ewcut) and $
;     (sdssispec.h_beta_sigma[0] lt sdss_sigmacut) and $
      (sdssispec.oii_3727[1] ne -2.0) and $
      (sdssispec.oiii_5007[1] ne -2.0) and $ ; nominal
      (sdssispec.nii_6584[1] ne -2.0) and $ ; nominal
      (sdssispec.h_alpha[1] ne -2.0)) ; nominal

    sel2 = where($
      (sdssispec.h_beta[0] gt sdss_hbcut) and $
      (sdssispec.h_beta[0]/sdssispec.h_beta[1] gt sdss_snrcut) and $ ; nominal cut
;     (sdssispec.h_beta_ew[0] gt sdss_ewcut) and $
;     (sdssispec.h_beta_sigma[0] lt sdss_sigmacut) and $
      (sdssispec.oii_3727[1] ne -2.0) and $
      (sdssispec.oiii_5007[1] ne -2.0) and $ ; nominal
      (sdssispec.nii_6584[1] ne -2.0) and $ ; nominal
      (sdssispec.h_alpha[1] ne -2.0) and $ ; nominal
      (sdssispec.oii_3727[0] gt 10^sdss_oiihbcut*sdssispec.h_beta[0]))
    rej2 = where($
      (sdssispec.h_beta[0] gt sdss_hbcut) and $
      (sdssispec.h_beta[0]/sdssispec.h_beta[1] gt sdss_snrcut) and $ ; nominal cut
;     (sdssispec.h_beta_ew[0] gt sdss_ewcut) and $
;     (sdssispec.h_beta_sigma[0] lt sdss_sigmacut) and $
      (sdssispec.oii_3727[1] ne -2.0) and $
      (sdssispec.oiii_5007[1] ne -2.0) and $ ; nominal
      (sdssispec.nii_6584[1] ne -2.0) and $ ; nominal
      (sdssispec.h_alpha[1] ne -2.0) and $ ; nominal
      (sdssispec.oii_3727[0] le 10^sdss_oiihbcut*sdssispec.h_beta[0]) and $
      (sdssispec.oii_3727[1] gt 0.0))
    help, sel1, sel2
    
; H-beta flux and EW vs redshift    
    psfile = qapath+'sdss_z_vs_hb.ps'
    im_plotconfig, 6, pos, psfile=psfile, xmargin=[1.7,0.4], $
      width=6.0, height=[3.0,3.0]

; luminosity    
    djs_plot, [0], [0], /nodata, position=pos[*,0], xsty=1, ysty=1, $
      xrange=zrange, yrange=scale*[1.05D38,9D41], xtitle='', xtickname=replicate(' ',10), $
      ytitle='L(H\beta) (10^{40} erg s^{-1})', /ylog
    djs_oplot, sdssispec[sel1].z, scale*sdss_dfactor[sel1]*sdssispec[sel1].h_beta[0], $
      psym=symcat(16), symsize=0.1, bin=20000L
    djs_oplot, sdssispec[rej1].z, scale*sdss_dfactor[rej1]*sdssispec[rej1].h_beta[0], $
      psym=symcat(9,thick=2), symsize=0.1, bin=5000L, color=im_color('dodger blue',101)
    djs_oplot, zaxis, scale*sdss_hbcut*4.0*!dpi*dluminosity(zaxis,/cm)^2, line=0, thick=6, color='red'
    im_legend, 'F(H\beta)>1\times10^{-16} erg s^{-1} cm^{-2}', $
      /right, /bottom, box=0, charsize=1.4, margin=0, line=0, thick=6

; EW
    djs_plot, [0], [0], /nodata, /noerase, position=pos[*,1], xsty=1, ysty=1, $
      xrange=zrange, yrange=[0.01,150], xtitle='Redshift', $
      ytitle='EW(H\beta) (\AA)', /ylog
    djs_oplot, sdssispec[sel1].z, sdssispec[sel1].h_beta_ew[0], $
      psym=symcat(16), symsize=0.1, bin=20000
    djs_oplot, sdssispec[rej1].z, sdssispec[rej1].h_beta_ew[0], $
      psym=symcat(9,thick=2), symsize=0.1, bin=5000L, color=im_color('dodger blue',101)
    im_plotconfig, /psclose

; [OII]/Hb and [OIII]/Hb vs redshift
    psfile = qapath+'sdss_z_vs_oiihb_oiiihb.ps'
    im_plotconfig, 6, pos, psfile=psfile, xmargin=[1.7,0.4], $
      width=6.0, height=[3.0,3.0], charsize=1.6

; [OII]/Hb    
    djs_plot, [0], [0], /nodata, position=pos[*,0], xsty=1, ysty=1, $
      xrange=zrange, yrange=[-1.5,1.0], xtitle='', xtickname=replicate(' ',10), $
      ytitle='log ([O II]/H\beta)'

    lim = where(sdssispec[sel1].oii_3727[1] eq -1.0,comp=det)
    djs_oplot, sdssispec[sel1[det]].z, alog10(sdssispec[sel1[det]].oii_3727[0]/sdssispec[sel1[det]].h_beta[0]), $
      psym=symcat(16), symsize=0.1, bin=20000L

    plotsym, 1.0, 0.1, thick=4
    djs_oplot, sdssispec[sel1[lim]].z, alog10(sdssispec[sel1[lim]].oii_3727_limit/sdssispec[sel1[lim]].h_beta[0]), $
      psym=8, color=im_color('orange red',101)
    djs_oplot, !x.crange, sdss_oiihbcut*[1,1], line=5, thick=6
    
; [OIII]/Hb
    djs_plot, [0], [0], /nodata, /noerase, position=pos[*,1], xsty=1, ysty=1, $
      xrange=zrange, yrange=[-1.8,1.5], xtitle='Redshift', $
      ytitle='log ([O III] \lambda5007/H\beta)'

    lim = where(sdssispec[sel2].oiii_5007[1] eq -1.0,comp=det)
    djs_oplot, sdssispec[sel2[det]].z, alog10(sdssispec[sel2[det]].oiii_5007[0]/$
      sdssispec[sel2[det]].h_beta[0]), psym=symcat(16), symsize=0.1, bin=20000L

    plotsym, 1.0, 0.1, thick=4
    djs_oplot, sdssispec[sel2[lim]].z, alog10(sdssispec[sel2[lim]].oiii_5007_limit/$
      sdssispec[sel2[lim]].h_beta[0]), psym=8, color=im_color('orange red',101), bin=5000L
    
    im_plotconfig, /psclose, psfile=psfile

; ------------------------------------------------------------
; U-V vs stellar mass for galaxies with and without upper limits
    agesispec = read_mz_sample(/mz_ispec)
    agesmass = read_mz_sample(/mz_mass)
    agesparent = read_mz_sample(/mz_ancillary)
    lim = where(agesispec.oiii_5007[1] eq -1.0,comp=det)

    mass = agesmass.mass_avg
    uv = agesparent.k_ubvrijhk_absmag_00[0]-agesparent.k_ubvrijhk_absmag_00[2]

    massrange = [8.0,12.0]
    uvrange = [0.4,2.3]
    
    psfile = qapath+'uv_vs_mass.ps'
    im_plotconfig, 0, pos, psfile=psfile, charsize=2.0, xmargin=[1.6,0.4], $
      height=6

    mzplot_scatterplot, mass[det], uv[det], position=pos, xsty=1, ysty=1, $
      xtitle=mzplot_masstitle(), ytitle='U - V', xrange=massrange, $
      yrange=uvrange;, npix=24
    djs_oplot, mass[lim], uv[lim], psym=6, sym=0.2, color=im_color('dodger blue',101)
    
    im_plotconfig, /psclose, psfile=psfile

; ------------------------------------------------------------
; EW(Hb), EW([OII]), EW([O III]) and EW(R23) vs redshift
    agesispec = read_mz_sample(/mzhii_ispec)
    
    psfile = qapath+'z_vs_ews.ps'
    im_plotconfig, 5, pos, psfile=psfile, height=3.0*[1,1], width=3.0*[1,1], $
      xmargin=[1.25,0.25], ymargin=[0.4,1.1], xspace=1.0, yspace=1.0, $
      charsize=1.5

    z = agesispec.z
    ewhb = agesispec.h_beta_ew[0]
    ewoii = agesispec.oii_3727_ew[0]
    ewoiii = agesispec.oiii_5007_ew[0]
    ewr23 = (ewoii+ocor*ewoiii)/ewhb

    xrange = [0.03,0.77]
    xtitle = 'Redshift'

; EW(Hb)    
    ytitle = 'EW(H\beta) (\AA)'
    yrange = [1.0,120.0]

    djs_plot, [0], [0], /nodata, position=pos[*,0], $
      xstyle=1, ystyle=1, xtitle=xtitle, ytitle=ytitle, $
      xrange=xrange, yrange=yrange, /ylog
    ages_oplot, z, ewhb
;   im_legend, '(a)', /left, /top, box=0, charsize=1.4, margin=0
    
; EW([OII])    
    ytitle = 'EW([O II]) (\AA)'
    yrange = [1.5,300.0]

    djs_plot, [0], [0], /nodata, /noerase, position=pos[*,1], $
      xstyle=1, ystyle=1, xtitle=xtitle, ytitle=ytitle, $
      xrange=xrange, yrange=yrange, /ylog
    ages_oplot, z, ewoii
;   im_legend, '(b)', /left, /top, box=0, charsize=1.4, margin=0
    
; EW([OIII])    
    ytitle = 'EW([O III] \lambda5007) (\AA)'
    yrange = [0.2,600.0]

    djs_plot, [0], [0], /nodata, /noerase, position=pos[*,2], $
      xstyle=1, ystyle=1, xtitle=xtitle, ytitle=ytitle, $
      xrange=xrange, yrange=yrange, /ylog
    ages_oplot, z, ewoiii
;   im_legend, '(c)', /left, /top, box=0, charsize=1.4, margin=0
    
; EW(R23)
    ytitle = 'EW(R_{23})'
;   ytitle = '[EW([O II])+EW([O III])]/EW(H\beta)'
    yrange = [0.8,25.0]

    djs_plot, [0], [0], /nodata, /noerase, position=pos[*,3], $
      xstyle=1, ystyle=1, xtitle=xtitle, ytitle=ytitle, $
      xrange=xrange, yrange=yrange, /ylog
    ages_oplot, z, ewr23
;   im_legend, '(d)', /left, /top, box=0, charsize=1.4, margin=0
    
    im_plotconfig, /psclose, psfile=psfile

return
end


;; ------------------------------------------------------------
;; redshift vs absolute magnitude and stellar mass, with the
;; completeness curves overlaid
;
;    psfile = pspath+'ages_redshift_vs_mr_mass'+suffix
;    im_plotconfig, 6, pos, psfile=psfile, height=3.0*[1,1], xmargin=[1.2,0.3], $
;      ymargin=[0.8,1.0], charsize=1.9, width=5.5
;
;    limits = mrdfits(mzpath+'mz_limits.fits.gz',1)
;
;    z = ageskcorr.z
;    mg = ageskcorr.ugriz_absmag[1]
;    gr = ageskcorr.ugriz_absmag[1]-ageskcorr.ugriz_absmag[2]
;    mass = ageskcorr.isedfit_mass
;    weight = ageskcorr.spec_weight
;    
;    xtitle = 'Redshift'
;    ytitle1 = textoidl('M_{0.1g} - 5 log (h_{70})')
;    ytitle2 = textoidl('log (M_{*}h_{70}^2/M'+sunsymbol()+')')
;
;    xrange = [0.03,0.77]
;    yrange1 = [-14.5,-24.0]
;    yrange2 = [7.5,12.0]
;
;; redshift vs ^{0.1} M_r
;    plot, [0], [0], /nodata, ysty=1, xsty=9, xrange=xrange, yrange=yrange1, $
;      xtitle='', ytitle=ytitle1, xtickname=replicate(' ',10), $
;      position=pos[*,0], yminor=4
;    axis, /xaxis, xsty=1, xrange=xrange, xtitle='Lookback Time (Gyr)', $
;      xtickv=getredshift(getage(0.0)-timelabel1), $
;      xticks=n_elements(timelabel1)-1L, xtickname=string(timelabel1,format='(I0)')
;    ages_oplot, z, mg
;;   djs_oplot, limits.zaxis, limits.mglim_50, line=0, thick=6.0, color='red'
;    djs_oplot, limits.zaxis, limits.mglim_75, line=5, thick=8.0, color='blue'
;    im_legend, ['75% !8K-correction!6 Completeness'], $
;      /right, /bottom, box=0, color='blue', line=5, pspacing=1.4, $
;      charsize=1.4, thick=8
;;   im_legend, ['50% K-correction Completeness','75% K-correction Completeness'], $
;;     /right, /bottom, box=0, color=['red','blue'], line=[0,5], pspacing=1.3, $
;;     charsize=1.5, thick=6.0
;; z vs stellar mass
;    plot, [0], [0], /nodata, /noerase, ysty=1, xsty=9, $
;      xrange=xrange, yrange=yrange2, xtitle=xtitle, ytitle=ytitle2, $
;      position=pos[*,1], yminor=5
;    ages_oplot, z, mass
;;   djs_oplot, limits.zaxis, limits.mmlim_50, line=0, thick=6, color='red'
;    djs_oplot, limits.zaxis, limits.mmlim_75, line=5, thick=8, color='blue'
;    im_legend, ['75% !8M/L!6 Completeness'], $
;      /right, /bottom, box=0, color='blue', line=5, pspacing=1.4, $
;      charsize=1.4, thick=8
;
;    im_plotconfig, /psclose
;


;;; ------------------------------------------------------------
;;; redshift vs absolute magnitude and stellar mass
;;    psfile = pspath+'z_vs_mb_mass'+suffix
;;    im_plotconfig, 6, pos, psfile=psfile, height=3.0*[1,1], xmargin=[1.3,0.2], $
;;      ymargin=[0.8,1.0], charsize=1.9, width=5.5
;;
;;; desired redshift grid
;;    sspzform = 5.0
;;    tauzform = 2.0
;;    csfzform = 2.0
;;    dz = 0.02 & minz = 0.0 & maxz = 1.0
;;    zz = (findgen((maxz-minz)/dz+1)*dz+minz)>0.01
;;    nz = n_elements(zz)
;;
;;; specify the models    
;;    modelspath = getenv('ISEDFIT_SFHGRID_DIR')+'/basemodels/bc03/'
;;    sspfits = 'chab_Z0.02_tau_00.0Gyr.fits.gz'
;;    taufits = 'chab_Z0.02_tau_03.0Gyr.fits.gz'
;;    csffits = 'chab_Z0.02_tau_100Gyr.fits.gz'
;;
;;    ssp = mrdfits(modelspath+sspfits,1,/silent)
;;    tau = mrdfits(modelspath+taufits,1,/silent)
;;    csf = mrdfits(modelspath+csffits,1,/silent)
;;    wave = csf.wave
;;
;;; interpolate the models at the appropriate ages, given ZFORM and the
;;; desired redshift grid
;;    sspflux = interpolate(ssp.flux,findex(ssp.age/1D9,getage(zz)-getage(sspzform)))
;;    tauflux = interpolate(tau.flux,findex(tau.age/1D9,getage(zz)-getage(tauzform)))
;;    csfflux = interpolate(csf.flux,findex(csf.age/1D9,getage(zz)-getage(csfzform)))
;;
;;; now given the apparent magnitude (I=20), compute the rest-frame
;;; luminosity and stellar mass at each redshift
;;    in_filterlist = 'ndwfs_I.par'
;;    out_filterlist = 'sdss_g0.par'
;;    sdss_band_shift = 0.1
;;    solarmag = k_solar_magnitudes(filterlist=out_filterlist,$
;;      band_shift=sdss_band_shift,/silent)
;;
;;    I20 = 19.95
;;    Ivega2ab = k_vega2ab(filterlist=in_filterlist,/kurucz,/silent)
;;    maggies = reform(10^(-0.4*(I20+Ivega2ab)),1,1)
;;
;;    mgvz = replicate({z: 0.0, ssp_Mg: 0.0, tau_Mg: 0.0, csf_Mg: 0.0, $
;;      ssp_mass: 0.0, tau_mass: 0.0, csf_mass: 0.0},nz)
;;    mgvz.z = zz
;;    for ii = 0L, nz-1L do begin
;;; SSP
;;       kk = im_simple_kcorrect(zz[ii],maggies,maggies*0.0+1.0,$
;;         in_filterlist,wave,sspflux[*,ii],band_shift=sdss_band_shift,$
;;         out_filterlist=out_filterlist,absmag=abs)
;;       mgvz[ii].ssp_Mg = abs
;;       mgvz[ii].ssp_mass = 10^(-0.4*(abs-solarmag))
;;; CSF
;;       kk = im_simple_kcorrect(zz[ii],maggies,maggies*0.0+1.0,$
;;         in_filterlist,wave,csfflux[*,ii],band_shift=sdss_band_shift,$
;;         out_filterlist=out_filterlist,absmag=abs)
;;       mgvz[ii].csf_Mg = abs
;;       mgvz[ii].csf_mass = 10^(-0.4*(abs-solarmag))
;;; TAU
;;       kk = im_simple_kcorrect(zz[ii],maggies,maggies*0.0+1.0,$
;;         in_filterlist,wave,tauflux[*,ii],band_shift=sdss_band_shift,$
;;         out_filterlist=out_filterlist,absmag=abs)
;;       mgvz[ii].tau_Mg = abs
;;       mgvz[ii].tau_mass = 10^(-0.4*(abs-solarmag))
;;    endfor
;;
;;; some other previously calculated limits    
;;    limits = mrdfits(mzpath+'mz_limits.fits.gz',1)
;;
;;; now make the plot    
;;    z = ageskcorr.z
;;    mg = ageskcorr.ugriz_absmag[1]
;;    gr = ageskcorr.ugriz_absmag[1]-ageskcorr.ugriz_absmag[2]
;;    mass = ageskcorr.isedfit_mass
;;    weight = ageskcorr.spec_weight
;;    
;;    xtitle = 'Redshift'
;;    ytitle1 = textoidl('M_{0.1g} - 5 log (h_{70})')
;;    ytitle2 = textoidl('log (M_{*}/h_{70}^{-2}M'+sunsymbol()+')') 
;;
;;    xrange = [0.03,0.77]
;;    yrange1 = [-15.3,-23.4]
;;    yrange2 = [7.7,11.9]
;;
;;; redshift vs ^{0.1} M_r
;;    plot, [0], [0], /nodata, ysty=1, xsty=1, xrange=xrange, yrange=yrange1, $
;;      xtitle='', ytitle=ytitle1, xtickname=replicate(' ',10), $
;;      position=pos[*,0], yminor=4
;;    ages_oplot, z, mg
;;    djs_oplot, mgvz.z, mgvz.csf_Mg, line=0, color=im_color('firebrick',101), thick=8
;;;   djs_oplot, mgvz.z, mgvz.tau_Mg, line=5, color=im_color('royal blue',101), thick=8
;;;   djs_oplot, mgvz.z, mgvz.ssp_Mg, line=5, color=im_color('royal blue',102), thick=8
;;;   djs_oplot, limits.zaxis, limits.mglim_50, line=0, thick=6.0, color='red'
;;;   djs_oplot, limits.zaxis, limits.mglim_75, line=5, thick=8.0, color='blue'
;;;   im_legend, ['75% !8K-correction!6 Completeness'], $
;;;     /right, /bottom, box=0, color='blue', line=5, pspacing=1.4, $
;;;     charsize=1.4, thick=8
;;;   im_legend, ['50% K-correction Completeness','75% K-correction Completeness'], $
;;;     /right, /bottom, box=0, color=['red','blue'], line=[0,5], pspacing=1.3, $
;;;     charsize=1.5, thick=6.0
;;; z vs stellar mass
;;    plot, [0], [0], /nodata, /noerase, ysty=1, xsty=9, $
;;      xrange=xrange, yrange=yrange2, xtitle=xtitle, ytitle=ytitle2, $
;;      position=pos[*,1], yminor=5
;;    ages_oplot, z, mass
;;    djs_oplot, mgvz.z, alog10(mgvz.csf_mass), line=0, $
;;      color=im_color('firebrick',101), thick=8
;;;   djs_oplot, mgvz.z, alog10(mgvz.tau_mass), line=5, $
;;;     color=im_color('royal blue',101), thick=8
;;;   djs_oplot, mgvz.z, alog10(mgvz.ssp_mass), line=5, $
;;;     color=im_color('royal blue',102), thick=8
;;;   djs_oplot, limits.zaxis, limits.mmlim_50, line=0, thick=6, color='red'
;;;   djs_oplot, limits.zaxis, limits.mmlim_75, line=5, thick=8, color='blue'
;;    im_legend, ['\psi(t)=const'], /right, /bottom, box=0, $
;;      pspacing=1.4, color='firebrick', line=0, charsize=1.5, thick=8
;;;   im_legend, ['75% !8M/L!6 Completeness'], $
;;;     /right, /bottom, box=0, color='blue', line=5, pspacing=1.4, $
;;;     charsize=1.4, thick=8
;;
;;    im_plotconfig, /psclose
;;


;;; ------------------------------------------------------------
;;; EW(Hb) vs S/N in all three lines showing the cuts used in
;;; BUILD_MZ_EMLINE_SAMPLE 
;;    psfile = pspath+'ewhb_vs_snr'+suffix
;;    im_plotconfig, 4, pos, psfile=psfile, yspace=0.0, $
;;      width=4.5, height=2.5*[1,1,1], xmargin=[1.6,0.4], $
;;      charsize=1.5
;;
;;    ewhb = allispec.h_beta_ew[0]
;;    snrhb = allispec.h_beta[0]/allispec.h_beta[1]
;;    snroii = allispec.oii_3727[0]/allispec.oii_3727[1]
;;    snroiii = allispec.oiii_5007[0]/allispec.oiii_5007[1]
;;
;;    ewhbcut = 1.0
;;    snrcut = 2.0
;;    keep1 = where((snrhb gt snrcut))
;;    keep2 = where((snrhb gt snrcut) and (snroii gt snrcut))
;;    keep3 = where((snrhb gt snrcut) and (snroii gt snrcut) and (snroiii gt snrcut))
;;;   im_plothist, allispec.z, bin=0.02, thick=6
;;;   im_plothist, allispec[keep1].z, bin=0.02, /over, color=im_color('orange')
;;;   im_plothist, allispec[keep2].z, bin=0.02, /over, color=im_color('red')
;;;   im_plothist, allispec[keep3].z, bin=0.02, /over, color=im_color('blue')
;;    
;;;   ww = (where((snrhb gt snrcut) and (snroii gt snrcut) and (snroiii lt snrcut)))[0:20]
;;;   ww = (where((snrhb gt snrcut) and (snroii lt snrcut) and (snroiii lt snrcut)))[0:20]
;;;   ww = where(snrhb lt 2 and ewhb gt 1.0 and allispec.isbroad eq 0)
;;;   ww = where(snrhb gt 5 and ewhb lt 1.0 and allispec.isbroad eq 0)
;;;   qaplot_ages_gandalf_specfit, allispec[ww], psfile='junk.ps', /solar
;;
;;    snrrange = [0.5,400.0]
;;    xrange = [0.2,150.0]
;;    xtitle = 'EW(H\beta) (\AA)'
;;
;;; EW(Hb)    
;;;   mzplot_scatterplot, ewhb, snrhb, position=pos[*,0], $
;;;     xstyle=1, ystyle=1, xtitle='', ytitle='S/N (H\beta)', $
;;;     xrange=xrange, yrange=snrrange, /ylog, /xlog, /ages, $
;;;     xtickname=replicate(' ',10)
;;    djs_plot, [0], [0], /nodata, position=pos[*,0], $
;;      xstyle=1, ystyle=1, xtitle='', ytitle='S/N (H\beta)', $
;;      xrange=xrange, yrange=snrrange, /ylog, /xlog, $
;;      xtickname=replicate(' ',10)
;;;   djs_oplot, ewhbcut*[1,1], 10^!y.crange
;;    djs_oplot, 10^!x.crange, snrcut*[1,1]
;;    ages_oplot, ewhb, snrhb[keep], symsize=0.25
;;;   ages_oplot, ewhb[keep], snrhb[keep], symsize=0.25
;;;   ages_oplot, ewhb[toss], snrhb[toss], symsize=0.4, color='red'
;;
;;; oii    
;;    djs_plot, [0], [0], /nodata, /noerase, position=pos[*,1], $
;;      xstyle=1, ystyle=1, xtitle='', ytitle='S/N ([OII] \lambda3727)', $
;;      xrange=xrange, yrange=snrrange, /ylog, /xlog, $
;;      xtickname=replicate(' ',10)
;;    ages_oplot, ewhb[keep], snroii[keep], symsize=0.25
;;    ages_oplot, ewhb[toss], snroii[toss], symsize=0.4, color='red'
;;;   djs_oplot, ewhbcut*[1,1], 10^!y.crange
;;
;;    plotsym, 1, 1.0
;;    oiilim = where(snroii lt 10^!y.crange[0],noiilim)
;;    djs_oplot, ewhb[oiilim], fltarr(noiilim)+1.0, psym=8, color='blue'
;;
;;; oiii    
;;    djs_plot, [0], [0], /nodata, /noerase, position=pos[*,2], $
;;      xstyle=1, ystyle=1, xtitle=xtitle, ytitle='S/N ([OIII] \lambda5007)', $
;;      xrange=xrange, yrange=snrrange, /ylog, /xlog
;;    ages_oplot, ewhb[keep], snroiii[keep], symsize=0.25
;;    ages_oplot, ewhb[toss], snroiii[toss], symsize=0.4, color='red'
;;;   djs_oplot, ewhbcut*[1,1], 10^!y.crange
;;
;;    plotsym, 1, 1.0
;;    oiiilim = where(snroiii lt 10^!y.crange[0],noiiilim)
;;    djs_oplot, ewhb[oiiilim], fltarr(noiiilim)+1.0, psym=8, color='blue'
;;    
;;    im_plotconfig, /psclose
;;



;;; ------------------------------------------------------------
;;; EW(Hb), EW([OII]), EW([O III]) and EW(R23) measured from our fluxed
;;; and unfluxed spectra
;;    psfile = pspath+'ages_ew_tests'+suffix
;;    im_plotconfig, 3, pos, psfile=psfile, xspace=[1.0,1.0], charsize=1.7
;;
;;    indx = where((allispec.h_beta_ew[0] gt 0.0) and (unfluxed.h_beta_ew[0] gt 0.0) and $
;;      (allispec.oii_3727_ew[0] gt 0.0) and (unfluxed.oii_3727_ew[0] gt 0.0) and $
;;      (allispec.oiii_5007_ew[0] gt 0.0) and (unfluxed.oiii_5007_ew[0] gt 0.0))
;;    
;;    z = allispec[indx].z
;;    ewha = allispec[indx].h_alpha_ew[0]
;;    ewhb = allispec[indx].h_beta_ew[0]
;;    ewoii = allispec[indx].oii_3727_ew[0]
;;    ewoiii = allispec[indx].oiii_5007_ew[0]
;;    ewr23 = (ewoii+ocor*ewoiii)/ewhb
;;
;;    ewoiicor = 1.2
;;    ewhbcor = 2.4
;;    ewhacor = 0.8
;;    ewoiiicor = 0.3
;;    
;;    uewha = unfluxed[indx].h_alpha_ew[0]
;;    uewhb = unfluxed[indx].h_beta_ew[0]
;;    uewoii = unfluxed[indx].oii_3727_ew[0]
;;    uewoiii = unfluxed[indx].oiii_5007_ew[0]
;;    uewr23 = ((uewoii+ewoiicor)+ocor*(uewoiii+ewoiiicor))/(uewhb+ewhbcor)
;;
;;; stats
;;    oiistats = im_stats(100*(ewoii-uewoii)/(0.5*(ewoii+uewoii)),sigrej=3.0)
;;    oiiistats = im_stats(100*(ewoiii-uewoiii)/(0.5*(ewoiii+uewoiii)),sigrej=3.0)
;;    hbstats = im_stats(100*(ewhb-uewhb)/(0.5*(ewhb+uewhb)),sigrej=3.0)
;;    hastats = im_stats(100*(ewha-uewha)/(0.5*(ewha+uewha)),sigrej=3.0)
;;    r23stats = im_stats(100*(ewr23-uewr23)/(0.5*(ewr23+uewr23)),sigrej=3.0)
;;;   oiistats = im_stats(ewoii-uewoii,sigrej=3.0)
;;;   oiiistats = im_stats(ewoiii-uewoiii,sigrej=3.0)
;;;   hbstats = im_stats(ewhb-uewhb,sigrej=3.0)
;;;   hastats = im_stats(ewha-uewha,sigrej=3.0)
;;;   r23stats = im_stats(ewr23-uewr23,sigrej=3.0)
;;
;;    splog, '[OII]:  '+$
;;      strtrim(string(oiistats.median_rej,format='(F12.4)'),2)+' ('+$
;;      strtrim(string(oiistats.mean_rej,format='(F12.4)'),2)+'+/-'+$
;;      strtrim(string(oiistats.sigma_rej,format='(F12.4)'),2)+')'
;;    splog, '[OIII]: '+$
;;      strtrim(string(oiiistats.median_rej,format='(F12.4)'),2)+' ('+$
;;      strtrim(string(oiiistats.mean_rej,format='(F12.4)'),2)+'+/-'+$
;;      strtrim(string(oiiistats.sigma_rej,format='(F12.4)'),2)+')'
;;    splog, 'Hb    : '+$
;;      strtrim(string(hbstats.median_rej,format='(F12.4)'),2)+' ('+$
;;      strtrim(string(hbstats.mean_rej,format='(F12.4)'),2)+'+/-'+$
;;      strtrim(string(hbstats.sigma_rej,format='(F12.4)'),2)+')'
;;    splog, 'Ha    : '+$
;;      strtrim(string(hastats.median_rej,format='(F12.4)'),2)+' ('+$
;;      strtrim(string(hastats.mean_rej,format='(F12.4)'),2)+'+/-'+$
;;      strtrim(string(hastats.sigma_rej,format='(F12.4)'),2)+')'
;;    splog, 'R23   : '+$
;;      strtrim(string(r23stats.median_rej,format='(F12.4)'),2)+' ('+$
;;      strtrim(string(r23stats.mean_rej,format='(F12.4)'),2)+'+/-'+$
;;      strtrim(string(r23stats.sigma_rej,format='(F12.4)'),2)+')'
;;    
;;; make the plot    
;;;   levels = errorf(0.5*(dindgen(2)+2)) ; 84.3%,96.7%
;;    
;;; EW(Hb)    
;;    xtitle = textoidl('EW(H\beta) (\AA) [Unfluxed]')
;;    ytitle = textoidl('EW(H\beta) (\AA) [Fluxed]')
;;    xrange = [0.2,99.0] & yrange = xrange
;;    xlogrange = alog10(xrange) & ylogrange = xlogrange
;;
;;    mzplot_scatterplot, /ages, alog10(uewhb), alog10(ewhb), position=pos[*,0], $
;;      xtitle=xtitle, ytitle=ytitle, xstyle=1, ystyle=1, xrange=xlogrange, yrange=ylogrange, $
;;      levels=levels;, ynpix=40, xnpix=40
;;    djs_oplot, logewaxis, logewaxis, line=0, thick=6.0
;;    djs_oplot, logewaxis, alog10(ewaxis+ewhbcor), line=5, color='dark red', thick=8.0
;;    legend, textoidl('2.4 \AA'), /right, /bottom, box=0, charsize=1.4, $
;;      margin=0, line=5, thick=8.0, color=im_color('dark red'), pspacing=1.5
;;
;;;   mzages_hogg_scatterplot, z, uewhb-ewhb, outpsym=symcat(agessym2), $
;;;     outsymsize=agespsize2, outcolor=im_color(agescolor2,1E8), $
;;;     xtitle='Redshift', ytitle=textoidl('\Delta'+'EW (\AA)'), $
;;;     charsize=charsize_1, xsty=1, ysty=1, xrange=xlogrange, yrange=ylogrange, $
;;;     /internal, ynpix=30, xnpix=30, /outliers, /noerase, $
;;;     position=[pos[2,0]-0.13,pos[1,0]+0.06,pos[2,0]-0.01,pos[1,0]+0.12]
;;    
;;; EW([OIII])    
;;    xtitle = textoidl('EW([O III]) (\AA) [Unfluxed]')
;;    ytitle = textoidl('EW([O III]) (\AA) [Fluxed]')
;;    xrange = [0.2,220.0] & yrange = xrange
;;    xlogrange = alog10(xrange) & ylogrange = xlogrange
;;
;;    mzplot_scatterplot, /ages, alog10(uewoiii), alog10(ewoiii), /noerase, position=pos[*,1], $
;;      xtitle=xtitle, ytitle=ytitle, xstyle=1, ystyle=1, xrange=xlogrange, yrange=ylogrange, $
;;      levels=levels;, ynpix=40, xnpix=40
;;    djs_oplot, logewaxis, logewaxis, line=0, thick=6.0
;;;   djs_oplot, logewaxis, alog10(ewaxis+ewoiiicor), line=5, color='dark red', thick=8.0
;;;   legend, textoidl('0.3 \AA'), /right, /bottom, box=0, charsize=1.4, $
;;;     margin=0, line=5, thick=8.0, color=im_color('dark red'), pspacing=1.5
;;
;;; EW([OII])    
;;    xtitle = textoidl('EW([O II]) (\AA) [Unfluxed]')
;;    ytitle = textoidl('EW([O II]) (\AA) [Fluxed]')
;;    xrange = [2.0,140.0] & yrange = xrange
;;    xlogrange = alog10(xrange) & ylogrange = xlogrange
;;
;;    mzplot_scatterplot, /ages, alog10(uewoii), alog10(ewoii), /noerase, position=pos[*,2], $
;;      xtitle=xtitle, ytitle=ytitle, xstyle=1, ystyle=1, xrange=xlogrange, yrange=ylogrange, $
;;      levels=levels;, ynpix=40, xnpix=40
;;    djs_oplot, logewaxis, logewaxis, line=0, thick=6.0
;;;   djs_oplot, logewaxis, alog10(ewaxis+ewoiicor), line=5, color='dark red', thick=linethick1
;;
;;; bottom panel: R23 residuals vs redshift    
;;    
;;;; EW(R23)
;;;    xtitle = textoidl('log EW(R_{23}) [Unfluxed]')
;;;    ytitle = textoidl('log EW(R_{23}) [Fluxed]')
;;;
;;;    xrange = [0.6,23.0] & yrange = xrange
;;;    xlogrange = alog10(xrange) & ylogrange = xlogrange
;;;
;;;    mzplot_scatterplot, /ages, alog10(uewr23), alog10(ewr23), /noerase, position=pos[*,3], $
;;;      xtitle=xtitle, ytitle=ytitle, xstyle=1, ystyle=1, xrange=xlogrange, yrange=ylogrange, $
;;;      levels=levels, ynpix=40, xnpix=40
;;;    djs_oplot, logewaxis, logewaxis, line=0, thick=6.0
;;;    im_legend, '(d)', /left, /top, box=0, charsize=1.4, margin=0
;;;;   im_legend, '(d) R_{23}', /left, /top, box=0, charsize=1.4, margin=0
;;    
;;    im_plotconfig, /psclose
;;
;;stop    
;;    

;;
;;; ------------------------------------------------------------
;;; EW(Hb) vs S/N in all three lines showing the cuts used in
;;; BUILD_MZ_EMLINE_SAMPLE 
;;    psfile = pspath+'ewhb_vs_snr'+suffix
;;    im_plotconfig, 4, pos, psfile=psfile, width=4.5, $
;;      height=2.5*[1,1,1], xmargin=[1.6,0.4], charsize=1.5, $
;;      yspace=0.8*[1.0,1.0]
;;
;;    ewhb = allispec.h_beta_ew[0]
;;    ewoii = allispec.oii_3727_ew[0]
;;    ewoiii = allispec.oiii_5007_ew[0]
;;
;;    snrhb = allispec.h_beta[0]/allispec.h_beta[1]
;;    snroii = allispec.oii_3727[0]/allispec.oii_3727[1]
;;    snroiii = allispec.oiii_5007[0]/allispec.oiii_5007[1]
;;
;;    hb = allispec.h_beta[0]
;;    oii = allispec.oii_3727[0]
;;    oiii = allispec.oiii_5007[0]
;;    
;;    ewcut = 0.0
;;    snrcut = 2.0
;;    keep1 = where((snrhb gt snrcut) and (ewhb gt ewcut))
;;    keep2 = where($
;;      (snrhb gt snrcut) and (ewhb gt ewcut) and $
;;      (snroii gt snrcut) and (ewoii gt ewcut))
;;    keep3 = where($
;;      (snrhb gt snrcut) and (ewhb gt ewcut) and $
;;      (snroii gt snrcut) and (ewoii gt ewcut) and $
;;      (snroiii gt snrcut) and (ewoiii gt ewcut))
;;    toss = where($
;;      (snrhb gt snrcut) and (ewhb gt ewcut) and $
;;      (snroii gt snrcut) and (ewoii gt ewcut) and $
;;      ((snroiii le snrcut) or (ewoiii le ewcut)))
;;
;;    check1 = where($
;;      ((snrhb lt snrcut) or (ewhb lt ewcut)) and $
;;      (snroii gt snrcut) and (ewoii gt ewcut) and $
;;      (snroiii gt snrcut) and (ewoiii gt ewcut))
;;    check2 = where($
;;      (snrhb gt snrcut) and (ewhb gt ewcut) and $
;;      ((snroii lt snrcut) or (ewoii lt ewcut)) and $
;;      (snroiii gt snrcut) and (ewoiii gt ewcut))
;;    
;;;   toss = where((snrhb gt snrcut) and (snroii gt snrcut) and (snroiii lt snrcut))
;;;   im_plothist, allispec.z, bin=0.02, thick=6
;;;   im_plothist, allispec[keep1].z, bin=0.02, /over, color=im_color('orange')
;;;   im_plothist, allispec[keep2].z, bin=0.02, /over, color=im_color('red')
;;;   im_plothist, allispec[keep3].z, bin=0.02, /over, color=im_color('blue')
;;    mv = agesparent.k_ubvrijhk_absmag_00[2]
;;    uv = agesparent.k_ubvrijhk_absmag_00[0]-agesparent.k_ubvrijhk_absmag_00[2]
;;
;;    djs_plot, mv, uv, psym=3, yr=[0,2.5], xr=[-26,-16]
;;;   djs_oplot, mv[keep3], uv[keep3], psym=6, sym=0.2, color='yellow'
;;    djs_oplot, mv[keep2], uv[keep2], psym=6, sym=0.2, color='cyan'
;;    djs_oplot, mv[toss], uv[toss], psym=6, sym=0.2, color='red'
;;
;;    djs_plot, ewhb[keep3], ewoii[keep3], psym=6, sym=0.2, xsty=3, ysty=3
;;    
;;    djs_plot, ewhb[keep3], ewoii[keep3], psym=6, sym=0.2, $
;;      xsty=3, ysty=3, /xlog, /ylog, xr=[0.1,100], yr=[0.1,100]
;;    djs_oplot, ewhb[toss], ewoii[toss], psym=6, sym=0.2, color='blue'
;;    ewaxis = range(0.1,100,500,/log)
;;    djs_oplot, ewaxis, 18*ewaxis-6
;;
;;; redshift vs mass
;;    djs_plot, agesparent[keep3].z, agesparent[keep3].k_mass, ysty=3, psym=6
;;    djs_oplot, agesparent[toss].z, agesparent[toss].k_mass, psym=6, color='cyan'
;;    
;;;   ww = where((uv[keep3] gt 2.3))
;;;   qaplot_ages_gandalf_specfit, allispec[keep3[ww]], psfile='junk.ps', /solar
;;    
;;;   ww = (where((snrhb gt snrcut) and (snroii gt snrcut) and (snroiii lt snrcut)))[0:20]
;;;   ww = (where((snrhb gt snrcut) and (snroii lt snrcut) and (snroiii lt snrcut)))[0:20]
;;;   ww = where(snrhb lt 2 and ewhb gt 1.0 and allispec.isbroad eq 0)
;;;   ww = where(snrhb gt 5 and ewhb lt 1.0 and allispec.isbroad eq 0)
;;;   qaplot_ages_gandalf_specfit, allispec[ww], psfile='junk.ps', /solar
;;
;;    snrrange = [0.5,400.0]
;;    xrange = [0.2,150.0]
;;
;;; EW(Hb)    
;;    djs_plot, [0], [0], /nodata, position=pos[*,0], $
;;      xstyle=1, ystyle=1, xtitle='EW(H\beta) (\AA)', $
;;      ytitle='S/N (H\beta)', xrange=xrange, yrange=snrrange, /ylog, /xlog
;;;   djs_oplot, ewhbcut*[1,1], 10^!y.crange
;;    djs_oplot, 10^!x.crange, snrcut*[1,1]
;;;   ages_oplot, ewhb, snrhb[keep], symsize=0.25
;;    ages_oplot, ewhb, snrhb, symsize=0.25
;;;   ages_oplot, ewhb[keep1], snrhb[keep1], symsize=0.25
;;;   ages_oplot, ewhb[toss], snrhb[toss], symsize=0.4, color='red'
;;
;;; oii    
;;    djs_plot, [0], [0], /nodata, /noerase, position=pos[*,1], $
;;      xstyle=1, ystyle=1, xtitle='EW([OII] \lambda3727) (\AA)', $
;;      ytitle='S/N ([OII] \lambda3727)', $
;;      xrange=xrange, yrange=snrrange, /ylog, /xlog
;;    ages_oplot, ewoii, snroii, symsize=0.25
;;    djs_oplot, 10^!x.crange, snrcut*[1,1]
;;;   ages_oplot, ewoii[keep1], snroii[keep1], symsize=0.25
;;;   ages_oplot, ewoii[toss], snroii[toss], symsize=0.4, color='red'
;;;   djs_oplot, ewhbcut*[1,1], 10^!y.crange
;;
;;;   plotsym, 1, 1.0
;;;   oiilim = where(snroii lt 10^!y.crange[0],noiilim)
;;;   djs_oplot, ewhb[oiilim], fltarr(noiilim)+1.0, psym=8, color='blue'
;;
;;; oiii    
;;    djs_plot, [0], [0], /nodata, /noerase, position=pos[*,2], $
;;      xstyle=1, ystyle=1, xtitle='EW([OIII] \lambda5007) (\AA)', $
;;      ytitle='S/N ([OIII] \lambda5007)', $
;;      xrange=xrange, yrange=snrrange, /ylog, /xlog
;;    ages_oplot, ewoiii, snroiii, symsize=0.25
;;;   djs_oplot, ewhbcut*[1,1], 10^!y.crange
;;
;;;   plotsym, 1, 1.0
;;;   oiiilim = where(snroiii[keep2] lt 10^!y.crange[0],noiiilim)
;;;   djs_oplot, ewhb[keep2[oiiilim]], fltarr(noiiilim)+1.0, psym=8, color='blue'
;;    
;;    im_plotconfig, /psclose
;;
;;stop    
;;    

;;; ------------------------------------------------------------
;;; distribution of global properties - SDSS
;;
;;    psfile = pspath+'sdss_histograms'+suffix
;;    im_plotconfig, 5, pos, psfile=psfile, height=[3.0,3.0], $
;;      width=[3.0,3.0], xmargin=[1.15,1.15], xspace=0.2, yspace=1.0, $
;;      charsize=1.7
;;
;;    zbinsize = 0.01
;;    absmagbinsize = 0.25
;;    colorbinsize = 0.04
;;    massbinsize = 0.15
;;    yfactor = 1.1
;;
;;    z_parent = sdssparent.z
;;    z = sdsskcorr.z    
;;
;;    mg_parent = sdssparent.ugriz_absmag[1]
;;    mg = sdsskcorr.ugriz_absmag[1]
;;
;;    gr_parent = sdssparent.ugriz_absmag[1]-sdssparent.ugriz_absmag[2]
;;    gr = sdsskcorr.ugriz_absmag[1]-sdsskcorr.ugriz_absmag[2]
;;
;;    mass_parent = sdssparent.mass
;;    mass = sdsskcorr.mass
;;    
;;    weight_parent = z_parent*0.0+1.0
;;    weight = z*0.0+1.0
;;
;;    splog, weighted_quantile(z_parent,weight_parent,quant=0.5), $
;;      im_weighted_mean(z_parent,weight_parent)
;;    splog, weighted_quantile(z,weight,quant=0.5), $
;;      im_weighted_mean(z,weight)
;;    
;;    histthick1 = 5.0
;;    fspacing = 0.04
;;
;;; ####################
;;; Redshift
;;; ####################
;;
;;    im_plothist, z_parent, bin=zbinsize, weight=weight_parent, $
;;      /noplot, xbin, ybin
;;
;;    xtitle = 'Redshift'
;;    ytitle = 'Number'
;;    xrange = [0.0,0.3]
;;    yrange = [0,max(ybin)*yfactor]
;;    
;;    djs_plot, [0], [0], /nodata, xsty=5, ysty=5, xrange=xrange, yrange=yrange, $
;;      position=pos[*,0], xtickname=replicate(' ',10), ytickname=replicate(' ',10)
;;    im_plothist, z_parent, bin=zbinsize, weight=weight_parent, $
;;      /overplot, /fill, fcolor=im_color('grey'), /fline, $
;;      forientation=45, fspacing=fspacing, color=im_color('grey')
;;    im_plothist, z, bin=zbinsize, weight=weight, /overplot, thick=histthick1
;;    djs_plot, [0], [0], /nodata, /noerase, xtitle=xtitle, ytitle=ytitle, $
;;      xsty=1, ysty=1, xrange=xrange, yrange=yrange, position=pos[*,0], $
;;      ytickinterval=1E4
;;;   im_legend, '(a)', /left, /top, box=0, margin=0
;;
;;; ####################
;;; ^{0.1} M_r
;;; ####################
;;    
;;    im_plothist, mg_parent, bin=absmagbinsize, xbin, ybin, $
;;      weight=weight_parent, /noplot
;;
;;    xtitle = textoidl('M_{0.1g} - 5 log (h_{70})')
;;    ytitle = 'Number'
;;    xrange = reverse([-16.1,-23.9])
;;    yrange = [0,max(ybin)*yfactor]
;;    
;;    djs_plot, [0], [0], /nodata, /noerase, xsty=5, ysty=5, xrange=xrange, yrange=yrange, $
;;      position=pos[*,1], xtickname=replicate(' ',10), ytickname=replicate(' ',10)
;;    im_plothist, mg_parent, bin=absmagbinsize, weight=weight_parent, $
;;      /overplot, /fill, fcolor=im_color('grey'), /fline, $
;;      forientation=45, fspacing=fspacing, color=im_color('grey')
;;    im_plothist, mg, bin=absmagbinsize, weight=weight, /overplot, $
;;      thick=histthick1
;;    djs_plot, [0], [0], /nodata, /noerase, xtitle=xtitle, ytitle='', xsty=1, $
;;      ysty=1, xrange=xrange, yrange=yrange, position=pos[*,1], $
;;      ytickname=replicate(' ',10)
;;    axis, /yaxis, yrange=yrange, ystyle=1, ytitle=textoidl(ytitle), $
;;      ytickinterval=1E4
;;;   im_legend, '(b)', /left, /top, box=0, margin=0
;;
;;; ####################
;;; ^{0.1} (g-r) color
;;; ####################
;;
;;    im_plothist, gr_parent, bin=colorbinsize, weight=weight_parent, $
;;      /noplot, xbin, ybin
;;
;;    xtitle = '^{0.1}(g - r)'
;;    ytitle = 'Number'
;;    xrange = [0.05,1.19]
;;    yrange = [0,max(ybin)*yfactor]
;;    
;;    djs_plot, [0], [0], /nodata, /noerase, xsty=5, ysty=5, xrange=xrange, yrange=yrange, $
;;      position=pos[*,2], xtickname=replicate(' ',10), ytickname=replicate(' ',10)
;;    im_plothist, gr_parent, bin=colorbinsize, weight=weight_parent, $
;;      /overplot, /fill, fcolor=im_color('grey'), /fline, forientation=45, $
;;      fspacing=fspacing, color=im_color('grey')
;;    im_plothist, gr, bin=colorbinsize, weight=weight, /overplot, $
;;      thick=histthick1
;;    djs_plot, [0], [0], /nodata, /noerase, xtitle=xtitle, xrange=xrange, ytitle=ytitle, $
;;      xsty=1, ysty=1, yrange=yrange, position=pos[*,2], ytickinterval=1E4
;;;   im_legend, '(c)', /left, /top, box=0, margin=0
;;
;;; ####################
;;; Stellar Mass
;;; ####################
;;
;;    im_plothist, mass_parent, bin=massbinsize, weight=weight_parent, $
;;      /noplot, xbin, ybin
;;
;;    xtitle = textoidl('log (M_{*}h_{70}^2/M'+sunsymbol()+')')
;;    ytitle = 'Number'
;;    xrange = [8.1,11.9]
;;    yrange = [0,max(ybin)*yfactor]
;;    
;;    djs_plot, [0], [0], /nodata, /noerase, xsty=5, ysty=5, xrange=xrange, yrange=yrange, $
;;      position=pos[*,3], xtickname=replicate(' ',10), ytickname=replicate(' ',10)
;;    im_plothist, mass_parent, bin=massbinsize, weight=weight_parent, $
;;      /overplot, /fill, fcolor=im_color('grey'), /fline, forientation=45, $
;;      fspacing=fspacing, color=im_color('grey')
;;    im_plothist, mass, bin=massbinsize, /overplot, weight=weight, thick=histthick1
;;    djs_plot, [0], [0], /nodata, /noerase, xtitle=xtitle, ytitle='', $
;;      xsty=1, ysty=1, xrange=xrange, yrange=yrange, $
;;      position=pos[*,3], yminor=3, ytickname=replicate(' ',10)
;;    axis, /yaxis, yrange=yrange, ystyle=1, ytitle=textoidl(ytitle), $
;;      ytickinterval=1E4
;;;   im_legend, '(d)', /left, /top, box=0, margin=0
;;
;;    im_plotconfig, /psclose
;;
