function get_rhosfr, data, rhosfr_err=rhosfr_err, local=local
    if keyword_set(local) then begin
       rhosfr = data.rhosfr_local
       rhosfr_err = data.rhosfr_local_err
    endif else begin
       rhosfr = data.rhosfr
       rhosfr_err = data.rhosfr_err
    endelse
; take the log        
    rhosfr_err = rhosfr_err/rhosfr/alog(10)
    rhosfr = alog10(rhosfr)
return, rhosfr
end

pro oplot_sfrd_model, zform=zform, zcut=zcut, log=log, keynote=keynote
; render a shaded region showing the minimum and maximum allowable
; SFHs, for the Salpeter IMF
    sfrd = get_sfrd_evolution(zaxis=zaxis,$
      zform=zform,zcut=zcut,dz=0.05)
    minsfrd = get_sfrd_evolution(zaxis=zaxis,$
      zform=zform,zcut=zcut,dz=0.05,/sfrdmin)
    maxsfrd = get_sfrd_evolution(zform=zform,$
      zcut=zcut,dz=0.05,/sfrdmax)
    if keyword_set(log) then zz = alog10(1+zaxis) else zz = zaxis
    if keyword_set(log) then begin
       col = 'powder blue' 
;      if keyword_set(keynote) then col = 'wheat' else col = 'powder blue' 
       polyfill, [zz,reverse(zz)], [minsfrd,reverse(maxsfrd)], $
         /data, /fill, color=fsc_color(col,100), noclip=0
;      polyfill, [zz,reverse(zz)], [minsfrd,reverse(maxsfrd)], $
;        /data, /line_fill, orientation=45, spacing=0.1, $
;        color=fsc_color(col,100), noclip=0
;      polyfill, [zz,reverse(zz)], [minsfrd,reverse(maxsfrd)], $
;        /data, /line_fill, orientation=135, spacing=0.1, $
;        color=fsc_color(col,100), noclip=0
    endif else begin
       col = 'wheat'
       polyfill, [zz,reverse(zz)], [minsfrd,reverse(maxsfrd)], $
         /data, /fill, color=fsc_color(col,100), noclip=0
    endelse
;   djs_oplot, zz, sfrd, color='red', thick=5
;   djs_oplot, zz, minsfrd, color='red', thick=5
;   djs_oplot, zz, maxsfrd, color='red', thick=5
return
end    

pro oplot_data, lf24, wiphu, errthick1=errthick1, $
  errthick2=errthick2, log=log, keycolor=keycolor, $
  makelegend=makelegend, keynote=keynote, col=col
    if keyword_set(keynote) then col = 'black' else col = 'black'
; local SFRD
    rhosfr = get_rhosfr(lf24[0],rhosfr_err=rhosfr_err,/local)
    plots, [0.0], rhosfr[0], symsize=3.8, psym=symcat(4,thick=errthick1), $
      color=fsc_color(col,103)
;   oploterror, [0.0], rhosfr[0], [0.0], rhosfr_err[0], $
;     symsize=3.8, psym=symcat(4,thick=errthick1), $
;     errthick=errthick2, color=fsc_color('black',103), $
;     errcolor=fsc_color('black',103)
;; Shupe+98
;    rhosfr = get_rhosfr(shupe,rhosfr_err=rhosfr_err)
;    oploterror, shupe.z, rhosfr[0], shupe.zerr, rhosfr_err[0], $
;      symsize=3.5, psym=symcat(4,thick=errthick1), $
;      errthick=errthick2, color=fsc_color('red',103), $
;      errcolor=fsc_color('red',103)
; Wiphu+10 at 0.13<z<0.6 (Salpeter IMF)
    indx = lindgen(4)+1
    rhosfr = get_rhosfr(wiphu[indx],rhosfr_err=rhosfr_err)
    if keyword_set(log) then begin
       zz = alog10(1+wiphu[indx].z)
       zzerr = wiphu[indx].zerr/(1+wiphu[indx].z)/alog(10)
    endif else begin
       zz = wiphu[indx].z
       zzerr = wiphu[indx].zerr
    endelse
    oploterror, zz, rhosfr[0,*], zzerr, rhosfr_err[0,*], $
      symsize=3.0, psym=symcat(6,thick=errthick1), $
      errthick=errthick2, color=fsc_color(col,101), $
      errcolor=fsc_color(col,101)
; Magnelli+09
    indx = lindgen(2)+5
    rhosfr = get_rhosfr(wiphu[indx],rhosfr_err=rhosfr_err)
    if keyword_set(log) then begin
       zz = alog10(1+wiphu[indx].z)
       zzerr = wiphu[indx].zerr/(1+wiphu[indx].z)/alog(10)
    endif else begin
       zz = wiphu[indx].z
       zzerr = wiphu[indx].zerr
    endelse
    oploterror, zz, rhosfr[0,*], zzerr, rhosfr_err[0,*], $
      symsize=3.0, psym=symcat(9,thick=errthick1), $
      errthick=errthick2, color=fsc_color(col,102), $
      errcolor=fsc_color(col,102)
    if keyword_set(makelegend) then begin
       im_legend, ['Shupe+98','Magnelli+09','Rujopakarn+10'], $
         /right, /bottom, psym=[4,9,6], symsize=[2.2,1.8,1.8], $
         box=0, charsize=1.3, symthick=errthick1, margin=0, $
         textcolor=keycolor, color=col
    endif
return
end

pro oplot_hopkins, xrange=xrange, log=log, keycolor=keycolor, $
  makelegend=makelegend, keynote=keynote
; overplot the data from Hopkins+04,06

    hh = rsex(getenv('PAPERSPATH')+'/literature/data/06hopkins.sex')
;   hh = hh[where(hh.z lt xrange[1])]
    uv = where(strmatch(hh.indicator,'*UV*',/fold))
    ha = where(strmatch(hh.indicator,'*Ha*',/fold) or $
      strmatch(hh.indicator,'*Hb*',/fold) or $
      strmatch(hh.indicator,'*OII*',/fold))
    ir = where(strmatch(hh.indicator,'*IR*',/fold))
    rad = where(strmatch(hh.indicator,'*RADIO*',/fold) or $
      strmatch(hh.indicator,'*xray*',/fold))

    if keyword_set(log) then begin
       zz = alog10(1+hh.z)
       zzerr_lo = hh.zerr_lo/(1+hh.z)/alog(10)
       zzerr_hi = hh.zerr_hi/(1+hh.z)/alog(10)
    endif else begin
       zz = hh.z
       zzerr_lo = hh.zerr_lo
       zzerr_hi = hh.zerr_hi
    endelse

    if keyword_set(keynote) then begin
       uvcolor = 'blue'
       hacolor = 'lime green'
       ircolor = 'firebrick'
       radcolor = 'purple'
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
    if keyword_set(keynote) then begin
       errthick1 = 4
    endif else begin
       errthick1 = 1
    endelse
    symsize1 = 1.3
    nohat = 1

; UV    
    oploterror, zz[uv], hh[uv].sfrd, zzerr_lo[uv], $
      hh[uv].sfrderr_lo, color=fsc_color(uvcolor,10), errcolor=fsc_color(uvcolor,10), $
      errthick=errthick1, /lobar, psym=symcat(uvpsym,thick=errthick1), symsize=symsize1, /nohat
    oploterror, zz[uv], hh[uv].sfrd, zzerr_hi[uv], $
      hh[uv].sfrderr_hi, color=fsc_color(uvcolor,10), errcolor=fsc_color(uvcolor,10), $
      errthick=errthick1, /hibar, psym=symcat(uvpsym,thick=errthick1), symsize=symsize1, /nohat
; Ha
    oploterror, zz[ha], hh[ha].sfrd, zzerr_lo[ha], $
      hh[ha].sfrderr_lo, color=fsc_color(hacolor,11), errcolor=fsc_color(hacolor,11), $
      errthick=errthick1, /lobar, psym=symcat(hapsym,thick=errthick1), symsize=symsize1, /nohat
    oploterror, zz[ha], hh[ha].sfrd, zzerr_hi[ha], $
      hh[ha].sfrderr_hi, color=fsc_color(hacolor,11), errcolor=fsc_color(hacolor,11), $
      errthick=errthick1, /hibar, psym=symcat(hapsym,thick=errthick1), symsize=symsize1, /nohat
; IR
    oploterror, zz[ir], hh[ir].sfrd, zzerr_lo[ir], $
      hh[ir].sfrderr_lo, color=fsc_color(ircolor,12), errcolor=fsc_color(ircolor,12), $
      errthick=errthick1, /lobar, psym=symcat(irpsym,thick=errthick1), symsize=symsize1, /nohat
    oploterror, zz[ir], hh[ir].sfrd, zzerr_hi[ir], $
      hh[ir].sfrderr_hi, color=fsc_color(ircolor,12), errcolor=fsc_color(ircolor,12), $
      errthick=errthick1, /hibar, psym=symcat(irpsym,thick=errthick1), symsize=symsize1, /nohat
; Radio
    oploterror, zz[rad], hh[rad].sfrd, zzerr_lo[rad], $
      hh[rad].sfrderr_lo, color=fsc_color(radcolor,13), errcolor=fsc_color(radcolor,13), $
      errthick=errthick1, /lobar, psym=symcat(radpsym,thick=errthick1), symsize=symsize1, /nohat
    oploterror, zz[rad], hh[rad].sfrd, zzerr_hi[rad], $
      hh[rad].sfrderr_hi, color=fsc_color(radcolor,14), errcolor=fsc_color(radcolor,14), $
      errthick=errthick1, /hibar, psym=symcat(radpsym,thick=errthick1), symsize=symsize1, /nohat

    if keyword_set(makelegend) then begin
       label = ['UV','IR','H\alpha,[OII]','Radio,Xray']
       im_legend, label, /right, /bottom, box=0, psym=[uvpsym,irpsym,hapsym,radpsym], $
         symsize=symsize1, color=[uvcolor,ircolor,hacolor,radcolor], $
         spacing=1.7, charsize=1.3, margin=0, textcolor=keycolor
    endif

return
end    

pro plot_cosmicimf_madau, ps=ps, keynote=keynote
; jm10mar25ucsd - build the Madau plot and overlay some model SFHs

    cosmicimfpath = ages_path(/projects)+'cosmicimf/'
    paperpath = ages_path(/papers)+'cosmicimf/'
    if keyword_set(keynote) then paperpath = $
      getenv('RESEARCHPATH')+'/meetings/10apr_florida/keynote/'
    if keyword_set(ps) then suffix = '.ps' else suffix = '.eps'

; read the results
    zbins = cosmicimf_zbins(nzbins,/lf24)
    lf24 = mrdfits(cosmicimfpath+'lf24_fit.fits.gz',1)
    wiphu = mrdfits(cosmicimfpath+'lf24_wiphu.fits.gz',1)
    shupe = mrdfits(cosmicimfpath+'lf24_shupe.fits.gz',1)
    
; ---------------------------------------------------------------------------
; 1-panel Madau plot - 0<z<5

    psfile = paperpath+'madau'+suffix
    im_plotconfig, 0, pos, psfile=psfile, xmargin=[1.3,0.4], $
      width=6.8, height=5, keynote=keynote
    
    if keyword_set(keynote) then keycolor = djs_icolor('white')
    errthick1 = 8
    errthick2 = 5
    xrange = [-0.05,0.8]
    yrange = [-2.2,-0.3]
    ytitle = cosmicimf_rhotitle(/sfr)

    djs_plot, [0], [0], /nodata, position=pos, xsty=1, ysty=1, $
      xtitle='log (1+z)', ytitle=ytitle, xrange=xrange, $
      yrange=yrange
    im_legend, 'Salpeter IMF (0.1-100 M_{\odot})', /left, /bottom, $
      box=0, margin=0, charsize=1.4, textcolor=keycolor
    oplot_sfrd_model, zform=5.0, zcut=1.0, /log, keynote=keynote
    oplot_hopkins, xrange=xrange1, /log, keycolor=keycolor, $
      /makelegend, keynote=keynote
    oplot_data, lf24, wiphu, errthick1=errthick1, $
      errthick2=errthick2, /log, keycolor=keycolor, $
      keynote=keynote, col=col
    im_legend, ['Shupe+98','Magnelli+09','Rujopakarn+10'], $
      /left, /top, psym=[4,9,6], symsize=[2.0,1.6,1.6], $
      box=0, charsize=1.3, symthick=errthick1, margin=0, $
      textcolor=keycolor, color=col

    im_plotconfig, /psclose
    if keyword_set(keynote) then begin
       spawn, 'ps2pdf13 '+psfile+' '+repstr(psfile,suffix,'.pdf'), /sh
       rmfile, psfile
    endif

; ---------------------------------------------------------------------------
; 2-panel Madau plot - 0<z<5

    psfile = paperpath+'madau_2panel'+suffix
    im_plotconfig, 6, pos, psfile=psfile, keynote=keynote
    
    errthick1 = 8
    errthick2 = 5
    xrange1 = [-0.06,1.4]
    xrange2 = [-0.05,0.8]
    yrange = [-2.2,-0.3]
    ytitle = cosmicimf_rhotitle(/sfr)

; #########################
; top panel    
    djs_plot, [0], [0], /nodata, position=pos[*,0], xsty=1, ysty=1, $
      xtitle='Redshift z', ytitle='', xrange=xrange1, yrange=yrange
    oplot_sfrd_model, zform=1.2, zcut=1.2, keynote=keynote
    oplot_hopkins, xrange=xrange1, keycolor=keycolor, keynote=keynote
    oplot_data, lf24, wiphu, errthick1=errthick1, $
      errthick2=errthick2, keycolor=keycolor, /makelegend, $
      keynote=keynote
    
; #########################
; bottom panel
    djs_plot, [0], [0], /nodata, /noerase, position=pos[*,1], xsty=1, ysty=1, $
      xtitle='log (1+z)', ytitle='', xrange=xrange2, yrange=yrange
    oplot_sfrd_model, zform=5.0, zcut=1.0, /log, keynote=keynote
    oplot_hopkins, xrange=xrange1, /log, keycolor=keycolor, $
      /makelegend, keynote=keynote
    oplot_data, lf24, wiphu, errthick1=errthick1, $
      errthick2=errthick2, /log, keycolor=keycolor, $
      keynote=keynote

; y-title
    xyouts, pos[0,0]-0.1, (pos[1,0]-pos[3,1])/2.0+pos[3,1], $
      textoidl(ytitle), align=0.5, orientation=90, /normal, $
      color=keycolor
    
    im_plotconfig, /psclose
    if keyword_set(keynote) then begin
       spawn, 'ps2pdf13 '+psfile+' '+repstr(psfile,suffix,'.pdf'), /sh
       rmfile, psfile
    endif

return
end
