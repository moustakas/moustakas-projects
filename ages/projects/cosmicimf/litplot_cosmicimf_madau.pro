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

pro litplot_cosmicimf_madau, ps=ps, keynote=keynote
; jm10apr05ucsd - show the Madau plot and the stellar mass density
; evolution based on data from the literature for my Florida talk

    talkpath = getenv('RESEARCHPATH')+'/meetings/10apr_florida/'
    if keyword_set(keynote) then talkpath = talkpath+'keynote/'

; Madau plot - 0<z<5
    if keyword_set(ps) then suffix = '.ps' else suffix = '.eps'

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
;   good = where(sfrd_rhostar gt 0)
;   djs_oplot, alog10(1+zaxis2), alog10(interpol(sfrd_rhostar[good],zaxis1[good],zaxis2)), $
;     line=0, color=keycolor, thick=8

;   label = ['\rho_{*}(z)','!MI!N\rho_{SFR}(z)']
;   im_legend, label, /right, /top, box=0, line=[5,0], $
;     thick=10, charsize=1.3, margin=0, pspacing=1.8, spacing=2.2, $
;     color=keycolor, textcolor=keycolor, charthick=3.0

    im_plotconfig, /psclose
    if keyword_set(keynote) then begin
       spawn, 'ps2pdf13 '+psfile+' '+repstr(psfile,suffix,'.pdf'), /sh
       rmfile, psfile
    endif

return
end
