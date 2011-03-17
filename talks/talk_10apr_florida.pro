function fit_hopkins, zaxis, sfrd_rhostar=sfrd_rhostar
; fit the Hopkins and Wilkins data
    hh = rsex(getenv('PAPERSPATH')+'/literature/data/06hopkins.sex')

    sfrd = 10.0^hh.sfrd
    sfrderr = total([[hh.sfrderr_lo],[hh.sfrderr_hi]],2)/2.0*$
      hh.sfrd*alog(10.0)*0.0+0.1 ; NO WEIGHTING!

    parinfo = replicate({value: 1.0, limited: [1,0], limits: [1E-8,0.0]},4)
    parinfo.value = [0.014,0.11,1.4,2.2]
    pp = mpfitexpr('(P[0]+P[1]*X)/(1.0+(X/P[2])^P[3])',$
      hh.z,sfrd,sfrderr,parinfo=parinfo,/quiet)
;   pp = [0.0166,0.1848,1.9474,2.6316] ; Cole et al. 2001
;   pp = [0.014,0.11,1.4,2.2]          ; Wilkins et al. 2008, for h=0.7
    splog, 'SFR density coefficients ', pp

    sfrdfit = (pp[0]+pp[1]*zaxis)/(1.0+(zaxis/pp[2])^pp[3])

; integrate the SFRD curve
    ageaxis = getage(0.0)-getage(zaxis)
    sfrd_rhostar = (1-0.28)*im_integral(ageaxis,sfrdfit*1D9,$ ; M_sun/Mpc^3
      ageaxis,replicate(max(ageaxis),n_elements(ageaxis)))

return, sfrdfit
end

function fit_wilkins, zaxis, sfrdfit=sfrdfit

    ww = rsex(getenv('PAPERSPATH')+'/literature/data/08wilkins.sex')
    zz = total([[ww.zmin],[ww.zmax]],2)/2.0

    parinfo = replicate({value: 1.0, limited: [1,0], limits: [1E-8,0.0]},3)
    parinfo.value = [0.0023,0.68,1.2]
;   cc = [0.0023,0.68,1.2] ; Wilkins et al. (2008)
    cc = mpfitexpr('(P[0]*exp(-P[1]*X^P[2]))',zz,$
      ww.omega,ww.omega_err*0.0+1.0,parinfo=parinfo,/quiet)
    splog, 'Stellar mass density coefficients ', cc

    rhostarfit = 1.36D11*cc[0]*exp(-cc[1]*zaxis^cc[2]) ; [M_sun Mpc^-3]

return, rhostarfit
end

pro oplot_wilkins, xrange=xrange, keycolor=keycolor
; overplot the data from Wilkins+08

    ww = rsex(getenv('PAPERSPATH')+'/literature/data/08wilkins.sex')

    rhostar = ww.omega*1.36D11 ; [M_sun Mpc^-3 for h=0.7]
    rhostarerr = ww.omega_err*1.36D11 ; [M_sun Mpc^-3 for h=0.7]

    rhostarerr = rhostarerr/rhostar/alog(10)
    rhostar = alog10(rhostar)
    
    zz = total([[ww.zmin],[ww.zmax]],2)/2.0
    zzerr = (ww.zmax-ww.zmin)/2.0
    zzerr = zzerr/(1+zz)/alog(10)
    zz = alog10(1+zz)
    
    wcolor = 'orange'
    wpsym = 16
    errthick1 = 4
    symsize1 = 1.3
    nohat = 1

    oploterror, zz, rhostar, zzerr, rhostarerr, $
      color=fsc_color(wcolor,10), errcolor=fsc_color(wcolor,10), $
      errthick=errthick1, psym=symcat(wpsym), symsize=symsize1, /nohat

return
end    

pro oplot_hopkins, xrange=xrange, keycolor=keycolor, $
  keynote=keynote
; overplot the data from Hopkins+04,06

    hh = rsex(getenv('PAPERSPATH')+'/literature/data/06hopkins.sex')
    uv = where(strmatch(hh.indicator,'*UV*',/fold))
    ha = where(strmatch(hh.indicator,'*Ha*',/fold) or $
      strmatch(hh.indicator,'*Hb*',/fold) or $
      strmatch(hh.indicator,'*OII*',/fold))
    ir = where(strmatch(hh.indicator,'*IR*',/fold))
    rad = where(strmatch(hh.indicator,'*RADIO*',/fold) or $
      strmatch(hh.indicator,'*xray*',/fold))

    zz = alog10(1+hh.z)
    zzerr_lo = hh.zerr_lo/(1+hh.z)/alog(10)
    zzerr_hi = hh.zerr_hi/(1+hh.z)/alog(10)

; make the plot    
    if keyword_set(keynote) then begin
       uvcolor = 'cyan'
       hacolor = 'khaki'
       ircolor = 'tomato'
       radcolor = 'orchid'
    endif else begin
       uvcolor = 'royal blue'
       hacolor = 'forest green'
       ircolor = 'red'
       radcolor = 'purple'
    endelse

    uvpsym = 9
    hapsym = 4
    irpsym = 15
    radpsym = 5
    errthick1 = 4
    symsize1 = 1.3
    nohat = 1

; UV    
    oploterror, zz[uv], hh[uv].sfrd, zzerr_lo[uv], $
      hh[uv].sfrderr_lo, color=fsc_color(uvcolor,10), errcolor=fsc_color(uvcolor,10), $
      errthick=errthick1, /lobar, psym=symcat(uvpsym), symsize=symsize1, /nohat
    oploterror, zz[uv], hh[uv].sfrd, zzerr_hi[uv], $
      hh[uv].sfrderr_hi, color=fsc_color(uvcolor,10), errcolor=fsc_color(uvcolor,10), $
      errthick=errthick1, /hibar, psym=symcat(uvpsym), symsize=symsize1, /nohat
; Ha
    oploterror, zz[ha], hh[ha].sfrd, zzerr_lo[ha], $
      hh[ha].sfrderr_lo, color=fsc_color(hacolor,11), errcolor=fsc_color(hacolor,11), $
      errthick=errthick1, /lobar, psym=symcat(hapsym), symsize=symsize1, /nohat
    oploterror, zz[ha], hh[ha].sfrd, zzerr_hi[ha], $
      hh[ha].sfrderr_hi, color=fsc_color(hacolor,11), errcolor=fsc_color(hacolor,11), $
      errthick=errthick1, /hibar, psym=symcat(hapsym), symsize=symsize1, /nohat
; IR
    oploterror, zz[ir], hh[ir].sfrd, zzerr_lo[ir], $
      hh[ir].sfrderr_lo, color=fsc_color(ircolor,12), errcolor=fsc_color(ircolor,12), $
      errthick=errthick1, /lobar, psym=symcat(irpsym), symsize=symsize1, /nohat
    oploterror, zz[ir], hh[ir].sfrd, zzerr_hi[ir], $
      hh[ir].sfrderr_hi, color=fsc_color(ircolor,12), errcolor=fsc_color(ircolor,12), $
      errthick=errthick1, /hibar, psym=symcat(irpsym), symsize=symsize1, /nohat
; Radio
    oploterror, zz[rad], hh[rad].sfrd, zzerr_lo[rad], $
      hh[rad].sfrderr_lo, color=fsc_color(radcolor,13), errcolor=fsc_color(radcolor,13), $
      errthick=errthick1, /lobar, psym=symcat(radpsym), symsize=symsize1, /nohat
    oploterror, zz[rad], hh[rad].sfrd, zzerr_hi[rad], $
      hh[rad].sfrderr_hi, color=fsc_color(radcolor,14), errcolor=fsc_color(radcolor,14), $
      errthick=errthick1, /hibar, psym=symcat(radpsym), symsize=symsize1, /nohat

    label = ['UV','IR','H\alpha,[OII]','Radio,Xray']
    im_legend, label, /right, /bottom, box=0, psym=[uvpsym,irpsym,hapsym,radpsym], $
      symsize=symsize1, color=[uvcolor,ircolor,hacolor,radcolor], $
      spacing=1.7, charsize=1.3, textcolor=keycolor, margin=0

return
end    

pro talk_10apr_florida, ps=ps, keynote=keynote
; jm10apr06ucsd - miscellaneous plots for my 2010/Apr Florida talk 

    common talkplots, primus1, kcorr1
    
    cosmicimfpath = ages_path(/projects)+'cosmicimf/'
    talkpath = getenv('RESEARCHPATH')+'/meetings/10apr_florida/'
    if keyword_set(keynote) then talkpath = talkpath+'keynote/'
    if keyword_set(ps) then suffix = '.ps' else suffix = '.eps'
    if keyword_set(keynote) then keycolor = djs_icolor('white')

    sample = read_cosmicimf_sample()

    if (n_elements(primus1) eq 0) or (n_elements(kcorr1) eq 0) then begin
       primus1 = mrdfits(talkpath+'primus_allz.fits.gz',1)
       kcorr1 = mrdfits(talkpath+'primus_allkcorr.fits.gz',1)
       spherematch, primus1.ra, primus1.dec, kcorr1.ra, $
         kcorr1.dec, 1.0/3600.0, m1, m2
       primus1 = primus1[m1]
       kcorr1 = kcorr1[m2]
    endif
    good = where(strtrim(primus1.zprimus_class,2) eq 'GALAXY' and $
      (kcorr1.k_mobs_sdss[3] lt 22.5))
;     (primus1.zprimus_zwarning eq 0))
    primus = primus1[good]
    kcorr = kcorr1[good]

; ---------------------------------------------------------------------------
; redshift histogram for PRIMUS
    psfile = talkpath+'primus_zhist'+suffix
    im_plotconfig, 0, pos, psfile=psfile, xmargin=[1.3,0.4], $
      width=6.8, height=4.5, keynote=keynote, charthick=4.0
    if keyword_set(keynote) then keycolor = djs_icolor('white') else $
      keycolor = djs_icolor('default')

    zbin = 0.01
    im_plothist, primus.zprimus, bin=zbin, xbin, ybin, /noplot

    xtitle = 'Redshift'
    ytitle = 'Number of Galaxies'
    xrange = [0.0,1.2]
    yrange = [0,max(ybin)*1.1]

; plot the redshift distribution for various absolute magnitude limits
    absmag = kcorr.k_absmag_sdss_05[2] ; M_r at band_shift=0.5
    mlim = [-21.5,-20.0,-19.0,-18.0,-17.0]
    if keyword_set(keynote) then $
      col = fsc_color(['cyan','green','red','grey','orange'],findgen(5)+50) else $
        col = fsc_color(['blue','green','red','grey','cyan'],findgen(5)+50)

    djs_plot, [0], [0], /nodata, xsty=5, ysty=5, xrange=xrange, yrange=yrange, $
      position=pos, xtickname=replicate(' ',10), ytickname=replicate(' ',10)
    these = where(absmag lt 0.0)
;   im_plothist, primus[these].zprimus, bin=zbin, /overplot, /fill, $
;     fcolor=djs_icolor('orange'), /fline, $
;     forientation=45, fspacing=fspacing, color=djs_icolor('orange')

    for ii = 0, n_elements(mlim)-1 do begin
       these = where((absmag gt mlim[ii]) and (absmag lt 0))
       im_plothist, primus[these].zprimus, bin=zbin, /overplot, /fill, $
         fcolor=djs_icolor(col[ii]), /fline, $
         forientation=45, fspacing=fspacing, color=djs_icolor(col[ii])
    endfor
    djs_plot, [0], [0], /nodata, /noerase, xtitle=xtitle, ytitle=ytitle, $
      xsty=1, ysty=1, xrange=xrange, yrange=yrange, position=pos
    im_legend, 'PRIMUS i_{AB}<22.5', /left, /top, box=0, textcolor=keycolor, $
      charsize=1.5, charthick=3.5, spacing=2.3, margin=0
    label = ['M_{0.5r}>'+string(mlim,format='(F5.1)')]
    im_legend, label, /right, /top, box=0, textcolor=col, $
      charsize=1.3, charthick=3.5, spacing=2.3, margin=0
    
    im_plotconfig, /psclose
    if keyword_set(keynote) then begin
       spawn, 'ps2pdf13 '+psfile+' '+repstr(psfile,'.eps','.pdf'), /sh
       rmfile, psfile
    endif

stop    
    
; ---------------------------------------------------------------------------
; show an example SED fit
    
    
    
; ---------------------------------------------------------------------------
; stellar mass vs redshift for PRIMUS
    psfile = talkpath+'primus_mass_vs_redshift'+suffix
    im_plotconfig, 0, pos, psfile=psfile, xmargin=[1.3,0.4], $
      width=6.8, height=4.5, keynote=keynote

    xtitle = 'Redshift'
    ytitle = cosmicimf_masstitle()
    xrange = [0.05,1.2]
    yrange = [8,12.5]

    mass = alog10(kcorr.k_kcorrect_mass) - alog10(0.7^2) + 0.26
;   mass = reform(k_sdss_bell(reform(kcorr.primus_absmag_sdss[0,*]))) ; stellar mass
;   mass = alog10(bmass) - alog10(0.7^2) + 0.14 ; h=1-->0.7; diet Salpeter-->Salpeter

    hogg_scatterplot, kcorr.primus_redshift, mass, position=pos, $
      xsty=1, ysty=1, xrange=xrange, yrange=yrange, /outliers, $
      /internal, outcolor=fsc_color('dodger blue',101), $
      xtitle=xtitle, ytitle=ytitle, levels=[0.5,0.75,0.9]
    
;   djs_plot, [0], [0], /nodata, xsty=1, ysty=1, xrange=xrange, $
;     yrange=yrange, position=pos
;   djs_oplot, kcorr.primus_redshift, mass, psym=3
    
    im_plotconfig, /psclose
    if keyword_set(keynote) then begin
       spawn, 'ps2pdf13 '+psfile+' '+repstr(psfile,'.eps','.pdf'), /sh
       rmfile, psfile
    endif

; ---------------------------------------------------------------------------
; pretty plot showing the Rieke+09 infrared SEDs    
    rr = read_09rieke()
    nrr = n_elements(rr.lir)
    loadct, 13
    col = range(1,255,nrr)
;   col = ['dark grey','steel blue','navy','cyan','green','orange',$
;     'firebrick',''
    
    psfile = talkpath+'rieke'+suffix
    im_plotconfig, 0, pos, psfile=psfile, keynote=keynote, $
      xmargin=[1.4,0.3], width=6.8, height=5

    djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
      xrange=[4,2E3], yrange=[2E-4,1.3], xtitle='Wavelength (\mu'+'m)', $
      ytitle='Relative flux (arbitrary units)', /xlog, /ylog
;   legend, 'Rieke+09', /left, /top, box=0, charsize=1.6
    wave = rr.wave/1D4 ; [micron]
    for ii = 0, nrr-1 do begin
       flux = rr.flux[*,ii]*rr.wave^2/im_light(/ang)
       if (ii eq 0) then norm = max(flux)
       scale = rr.lir[0]/rr.lir[ii]/norm
       djs_oplot, wave, scale*flux, color=col[ii]
    endfor
;   if keyword_set(keynote) then col[0] = djs_icolor('white')
    
    label = strtrim(string(alog10(rr.lir),format='(F12.2)'),2)
    legend, label[0:6], /right, /bottom, box=0, $
      charsize=1.3, color=col[0:6], textcolor=col[0:6], $
      line=0, pspacing=1.5, position=[100,3E-4], /data
    legend, label[7:nrr-1], /right, /bottom, box=0, $
      charsize=1.3, color=col[7:nrr-1], textcolor=col[7:nrr-1], $
      line=0, pspacing=1.5, position=[600,3E-4], /data

    im_plotconfig, /psclose
    if keyword_set(keynote) then begin
       spawn, 'ps2pdf13 '+psfile+' '+repstr(psfile,'.eps','.pdf'), /sh
       rmfile, psfile
    endif

; ---------------------------------------------------------------------------
; redshift histogram for AGES
    psfile = talkpath+'ages_zhist'+suffix
    im_plotconfig, 0, pos, psfile=psfile, xmargin=[1.3,0.4], $
      width=6.8, height=4.5, keynote=keynote

    zbin = 0.02
    im_plothist, sample.z, bin=zbin, xbin, ybin, $
      weight=sample.final_weight, /noplot

    xtitle = 'Redshift'
    ytitle = 'Number of Galaxies'
    xrange = [0.0,0.8]
    yrange = [0,max(ybin)*1.1]

    if keyword_set(keynote) then $
      histcolor = fsc_color('wheat',101) else $
      histcolor = djs_icolor('default')
    
    djs_plot, [0], [0], /nodata, xsty=5, ysty=5, xrange=xrange, yrange=yrange, $
      position=pos, xtickname=replicate(' ',10), ytickname=replicate(' ',10)
    im_plothist, sample.z, bin=zbin, weight=sample.final_weight, $
      /overplot, /fill, fcolor=djs_icolor('orange'), /fline, $
      forientation=45, fspacing=fspacing, color=djs_icolor('orange')
;   im_plothist, sample.z, bin=zbin, weight=sample.final_weight, $
;     /overplot, /fill, fcolor=djs_icolor('orange'), /fline, $
;     forientation=135, fspacing=fspacing, color=djs_icolor('orange')
    mips = where(sample.phot_mips24 gt 0.27)
    im_plothist, sample[mips].z, bin=zbin, weight=sample[mips].final_weight, $
      /overplot, thick=8.0, color=histcolor

    djs_plot, [0], [0], /nodata, /noerase, xtitle=xtitle, ytitle=ytitle, $
      xsty=1, ysty=1, xrange=xrange, yrange=yrange, position=pos
    legend, textoidl(['AGES I_{Vega}<20','MIPS f_{24}>0.27 mJy']), $
      /right, /top, box=0, textcolor=[djs_icolor('orange'),histcolor], $
      charsize=1.5, charthick=3.5, spacing=2.3, margin=0
    
    im_plotconfig, /psclose
    if keyword_set(keynote) then begin
       spawn, 'ps2pdf13 '+psfile+' '+repstr(psfile,'.eps','.pdf'), /sh
       rmfile, psfile
    endif

; ---------------------------------------------------------------------------
; redshift vs volume of various surveys
    surveys = ['PRIMUS','AGES','DEEP2','SDSS','TKRS','CNOC2',$
      'zCOSMOS-!cBright','zCOSMOS-!cDeep','COMBO17','2dFGRS']
    nsurvey = n_elements(surveys)
    info = {survey: '', nz: 0.0D, area: 0.0D, z12: [0.0,0.0], $
      vol: 0.0D, xs: 1.0, ys: 1.0}
    info = replicate(info,nsurvey)

    info.survey = surveys
    info.xs = [1.0,2.9,1.0,1.0,1.0,1.0,0.35,1.0,0.7,1.4]
    info.ys = [1.7,0.95,1.5,1.5,2.0,0.35,0.6,0.5,1.4,0.5]
    info.nz = [1.4E5,1.4E4,4E4,7.0E5,1440.0,2000.0,2E4,1E4,2.5E4,2.5E5]
    info.area = [10.0,7.8,3.0,7966.0,(10.0*16.0)/3600.0,1400.0/3600.0,1.7,1.0,0.78,2000.0]
    info.z12 = [[0.2,1.0],[0.01,0.8],[0.7,1.5],[0.01,0.25],[0.0,1.5],$
      [0.12,0.55],[0.1,1.2],[1.4,3.0],[0.2,1.1],[0.0,0.22]]
    for ii = 0, nsurvey-1 do info[ii].vol = $
      jhnvol(info[ii].z12[0],info[ii].z12[1])*info[ii].area*3600.0
    struct_print, info

    if keyword_set(keynote) then ccolor = 'cyan' else ccolor = 'dark red'

    xrange = [5E4,5E9]
    yrange = [500,1.8E6]

; focus on AGES
    info_ages = info[1:nsurvey-1]
    nsurvey_ages = nsurvey-1
    psfile = talkpath+'z_vs_volume_ages'+suffix
    im_plotconfig, 0, pos, xmargin=[1.3,0.2], charsize=2.2, $
      psfile=psfile, keynote=keynote
    djs_plot, [0], [0], /nodata, xrange=xrange, yrange=yrange, $
      /xlog, /ylog, /xsty, /ysty, xtitle='Comoving Volume (h_{70}^{-3} Mpc^{3})', $
      ytitle='Number of Redshifts', position=pos
    plots, info_ages[0].vol, info_ages[0].nz, psym=symcat(15,thick=5.0), $
      symsize=4, color=djs_icolor(ccolor)
    djs_oplot, info_ages[1:nsurvey_ages-1].vol, info_ages[1:nsurvey_ages-1].nz, $
      psym=symcat(4,thick=5.0), symsize=4.0
    xyouts, info_ages[0].vol*info_ages[0].xs, info_ages[0].nz*info_ages[0].ys, $
      info_ages[0].survey, charsize=1.6, align=0.5, color=djs_icolor(ccolor)
    for ii = 1, nsurvey_ages-1 do xyouts, info_ages[ii].vol*info_ages[ii].xs, $
      info_ages[ii].nz*info_ages[ii].ys, info_ages[ii].survey, charsize=1.6, align=0.5
    im_plotconfig, /psclose
    if keyword_set(keynote) then begin
       spawn, 'ps2pdf13 '+psfile+' '+repstr(psfile,'.eps','.pdf'), /sh
       rmfile, psfile
    endif

; focus on PRIMUS    
    psfile = talkpath+'z_vs_volume_primus'+suffix
    im_plotconfig, 0, pos, xmargin=[1.3,0.2], charsize=2.2, $
      psfile=psfile, keynote=keynote
    djs_plot, [0], [0], /nodata, xrange=xrange, yrange=yrange, $
      /xlog, /ylog, /xsty, /ysty, xtitle='Comoving Volume (h_{70}^{-3} Mpc^{3})', $
      ytitle='Number of Redshifts', position=pos
    djs_oplot, info[1:nsurvey-1].vol, info[1:nsurvey-1].nz, $
      psym=symcat(4,thick=5.0), symsize=4.0
    plots, info[0].vol, info[0].nz, psym=symcat(15,thick=5.0), $
      symsize=4, color=djs_icolor(ccolor)
    xyouts, info[0].vol*info[0].xs, info[0].nz*info[0].ys, info[0].survey, $
      charsize=1.6, align=0.5, color=djs_icolor(ccolor)
    for ii = 1, nsurvey-1 do xyouts, info[ii].vol*info[ii].xs, info[ii].nz*info[ii].ys, $
      info[ii].survey, charsize=1.6, align=0.5
    im_plotconfig, /psclose
    if keyword_set(keynote) then begin
       spawn, 'ps2pdf13 '+psfile+' '+repstr(psfile,'.eps','.pdf'), /sh
       rmfile, psfile
    endif

; ---------------------------------------------------------------------------
; Madau plot and the stellar mass density evolution based on data from
; the literature

; Madau plot - 0<z<5
    psfile = talkpath+'lit_madau_rhostar'+suffix
    im_plotconfig, 6, pos, psfile=psfile, keynote=keynote, $
      xmargin=[1.3,0.4], width=6.8, height=[3.4,3.4], $
      charsize=1.7
    if keyword_set(keynote) then keycolor = djs_icolor('white')

    zaxis1 = im_array(0.01,5.1,0.01)
    zaxis2 = im_array(0.01,4.5,0.01)

    xrange1 = [-0.03,0.8]
    yrange1 = [-2.2,-0.3]
    yrange2 = [7.0,9.3]

    ytitle1 = cosmicimf_rhotitle(/sfr)
    ytitle2 = cosmicimf_rhotitle()

; #########################
; top panel - madau plot
    djs_plot, [0], [0], /nodata, position=pos[*,0], xsty=1, ysty=1, $
      xtitle='', ytitle=ytitle1, xrange=xrange1, yrange=yrange1, $
      xtickname=replicate(' ',10)
    oplot_hopkins, xrange=xrange1, keycolor=keycolor, keynote=keynote
    sfrdfit = fit_hopkins(zaxis1,sfrd_rhostar=sfrd_rhostar)
    djs_oplot, alog10(1+zaxis1), alog10(sfrdfit), line=0, $
      color=keycolor, thick=8
    
; #########################
; bottom panel - mass density evolution
    djs_plot, [0], [0], /nodata, /noerase, position=pos[*,1], xsty=1, ysty=1, $
      xtitle='log (1+z)', ytitle=ytitle2, xrange=xrange1, yrange=yrange2
    oplot_wilkins, xrange=xrange1, keycolor=keycolor
    rhostarfit = fit_wilkins(zaxis2)
    djs_oplot, alog10(1+zaxis2), alog10(rhostarfit), line=5, color=keycolor, thick=8
    good = where(sfrd_rhostar gt 0)
    djs_oplot, alog10(1+zaxis2), alog10(interpol(sfrd_rhostar[good],zaxis1[good],zaxis2)), $
      line=0, color=keycolor, thick=8
    im_legend, 'Salpeter IMF (0.1-100 M_{\odot})', /left, /bottom, $
      box=0, margin=0, charsize=1.4, textcolor=keycolor

    label = ['\rho_{*}(z)','!MI!N\rho_{SFR}(z)']
    im_legend, label, /right, /top, box=0, line=[5,0], $
      thick=10, charsize=1.3, margin=0, pspacing=1.8, spacing=2.2, $
      color=keycolor, textcolor=keycolor, charthick=3.0

    im_plotconfig, /psclose
    if keyword_set(keynote) then begin
       spawn, 'ps2pdf13 '+psfile+' '+repstr(psfile,suffix,'.pdf'), /sh
       rmfile, psfile
    endif

return
end
