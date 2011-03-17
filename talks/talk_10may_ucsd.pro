pro ucsd_oplot_imf, mass, gamma=gamma, mcut=mcut, _extra=extra
    timf = imf(mass,gamma,mcut)
;   timf = 1E3*timf/total(timf,/double)
    timf = timf/interpol(timf,mass,1.0)
    good = where(timf ne 0)
    djs_oplot, alog10(mass[good]), alog10(timf[good]), $
      _extra=extra
return
end    

pro ucsd_ages_oplot, x, y, _extra=extra, keynote=keynote
    if keyword_set(keynote) then cc = 'cyan' else cc = 'grey'
    djs_oplot, x, y, psym=symcat(16), symsize=0.4, $
      color=cc, _extra=extra
return
end

pro ucsd_mzsdss_hogg_scatterplot, x, y, _extra=extra, keynote=keynote
; mass-metallicity plot wrapper on HOGG_SCATTERPLOT
    if keyword_set(keynote) then begin
       loadct, 1, /silent
       acolor = 'white'
       ccolor = 'black'
       outcolor = 'black'
       outsymsize = 0.15
    endif else begin
       loadct, 0, /silent
       acolor = 'black'
       ccolor = 'black'
       outcolor = 'grey'
       outsymsize = 0.3
    endelse
    levels = [0.01,0.05,0.1,0.3,0.6,0.9,0.95,0.99]
;   levels = [0.05,0.10,0.25,0.5,0.75,0.90,0.95]
    npix = 40
    eextra = struct_trimtags(extra,except=['POSITION'])
    pos = extra.position
    hogg_scatterplot, x, y, /outliers, outpsym=symcat(16), $
      outsymsize=outsymsize, outcolor=outcolor, _extra=eextra, $
      color=djs_icolor(acolor), ccolor=djs_icolor(ccolor), $
      position=pos, levels=levels, xnpix=npix, ynpix=npix, $
      exp=0.6, /internal
; overplot black axes
    eextra = struct_trimtags(extra,except=['XTITLE','YTITLE'])
    plot, [0], [0], /noerase, xsty=1, ysty=1, _extra=eextra, $
      color=djs_icolor(ccolor), ytitle='', xtitle='', $
      xtickname=replicate(' ',10), ytickname=replicate(' ',10)
    if keyword_set(keynote) then loadct, 0, /silent
return
end

pro talk_10may_ucsd, keynote=keynote
; jm10may10ucsd - miscellaneous plots for my 2010/May UCSD talk 

    mzpath = ages_path(/projects)+'mz/'
    talkpath = getenv('RESEARCHPATH')+'/talks/2010/10may_ucsd/'
    if keyword_set(keynote) then talkpath = talkpath+'keynote/'

    agesparent = read_mz_sample(/parent)
    agesmz = read_mz_sample(/mzhii_ancillary)

    sdssmz = read_mz_sample(/sdss,/mzhii_ancillary)
    sdssohnodust = read_mz_sample(/sdss,/mzhii_log12oh,/nodust)

; some basic plotting parameters    
    massrange1 = [8.1,12.0]
    massrange2 = [8.35,11.75]
    mgrange1 = [-16.1,-23.9]
    ohrange1 = [8.35,9.37]
    ohrange2 = [8.5,9.3]

    zaxis1 = im_array(0.0,0.8,0.01)

; ---------------------------------------------------------------------------    
; MZ relation from SDSS/cor
    psfile = talkpath+'ucsd_mzsdss.ps'
    im_plotconfig, 0, pos, psfile=psfile, charsize=2.0, $
      xmargin=[1.2,0.3], thick=4.0, keynote=keynote
    
    ucsd_mzsdss_hogg_scatterplot, sdssmz.mass_median, sdssmz.oh_median+randomn(seed,n_elements(sdssmz))*0.05, position=pos, $
;   s1 = mzlz_grab_info(sdssohnodust,sdssmz,/flux)
;   ucsd_mzsdss_hogg_scatterplot, s1.mass, s1.oh, position=pos, $
      xstyle=1, ystyle=1, xtitle=textoidl('log (M_{*}/M'+sunsymbol()+')'), ytitle=mz_ohtitle(), $
      xrange=[8.2,12.0], yrange=[8.0,9.5], keynote=keynote, $
      /noredraw

    im_plotconfig, /psclose
    if keyword_set(keynote) then begin
       spawn, 'ps2pdf13 '+psfile+' '+repstr(psfile,'.ps','.pdf'), /sh
       rmfile, psfile
    endif

stop    
    
; --------------------------------------------------
; plot some example IMFs
    psfile = talkpath+'salpeter_imf.ps'
    im_plotconfig, 0, pos, psfile=psfile, keynote=keynote, $
      xmargin=[1.3,0.4], width=[6.8], height=5.0

    if keyword_set(keynote) then keycolor = djs_icolor('white') else $
      keycolor = djs_icolor('default')

    mass1 = range(0.1D,100D,500,/log)
    mass2 = range(0.1D,120D,500,/log)

    ytitle = textoidl('log (\xi_{log m}) (dex^{-1})')
    xtitle = textoidl('log (m/M_{\odot})')
    
    plot, [0], [0], /nodata, position=pos, xsty=3, ysty=1, $
      xrange=alog10([0.1,100]), yrange=[-3,1.5], xtitle=xtitle, $
      ytitle=ytitle, charsize=2.2, color=keycolor
; xi = Number per unit logarithmic mass interval

; Salpeter
    salpcolor = 'red'
    ucsd_oplot_imf, mass1, gamma=-1.35, mcut=[0.1,100], line=0, $
      color=fsc_color(salpcolor,101), thick=10
;    gamma = [0.5,1.5,2.5]
;    line = [3,4,5]
;    if keyword_set(keynote) then $
;      color = ['cyan','wheat','orange'] else $
;        color = ['blue','dark green','orange']
;    for ii = 0, n_elements(gamma)-1 do ucsd_oplot_imf, mass2, gamma=[0.5,-gamma[ii]], $
;      mcut=[0.1,0.5,120], line=line[ii], color=fsc_color(color[ii],100), thick=10
;   im_legend, ['Salpeter','\Gamma='+string(gamma,format='(F3.1)')], $
;     line=[0,line], color=[salpcolor,color], thick=10, /left, $
;     /bottom, box=0, charsize=1.6, pspacing=1.8, textcolor=keycolor
;   im_legend, ['\xi_{logm}\propto'+'m^{-\Gamma}','Salpeter (\Gamma=1.35)'], $
;     /left, /bottom, box=0, line=[-1,0], charsize=2.0, pspacing=1.8, $
;     textcolor=keycolor
    im_legend, ['Salpeter (\Gamma=1.35)'], /left, /bottom, box=0, $
      line=0, color=salpcolor, thick=10, charsize=2.0, pspacing=1.8, $
      textcolor=keycolor
;   im_legend, ['Salpeter','\Gamma='+string(gamma,format='(F3.1)')], $
;     line=[0,line], color=[salpcolor,color], thick=10, /left, $
;     /bottom, box=0, charsize=1.6, pspacing=1.8, textcolor=keycolor

    im_plotconfig, /psclose
    if keyword_set(keynote) then begin
       spawn, 'ps2pdf13 '+psfile+' '+repstr(psfile,'.ps','.pdf'), /sh
       rmfile, psfile
    endif

    
stop    
    
; ------------------------------------------------------------
; redshift vs absolute magnitude and stellar mass

    psfile = talkpath+'ages_redshift_vs_mr_mass.ps'
    im_plotconfig, 6, pos, psfile=psfile, height=3.0*[1,1], xmargin=[1.3,0.2], $
      ymargin=[0.8,1.0], charsize=1.8, width=5.5, keynote=keynote

    if keyword_set(keynote) then keycolor = djs_icolor('white') else $
      keycolor = djs_icolor('default')

    zbin = 0.02
    im_plothist, agesparent.z, bin=zbin, xbin, ybin, $
      weight=agesparent.final_weight, /noplot
    hrange = [0,max(ybin)*1.1]

    htitle = 'Number'
    xtitle = 'Redshift'
    ytitle1 = mz_mgtitle()
    ytitle2 = mz_masstitle()

    xrange = [0.03,0.77]
    yrange1 = [-15.0,-23.2]
    yrange2 = [7.9,11.7]

; redshift vs M_{0.1g}
    plot, [0], [0], /nodata, ysty=1, xsty=1, xrange=xrange, yrange=yrange1, $
      xtitle='', ytitle=ytitle1, xtickname=replicate(' ',10), $
      position=pos[*,0], yminor=4, color=keycolor
    ucsd_ages_oplot, agesmz.z, agesmz.k_ugriz_absmag_01[1], keynote=keynote
;   djs_oplot, agesparent.z, agesparent.k_ugriz_absmag_01[1], ps=3

;; inset redshift histogram
;    hpos = [0.55,0.58,0.95,0.72]
;    plot, [0], [0], /nodata, /noerase, xsty=5, ysty=5, xrange=xrange, $
;      yrange=hrange, position=hpos, color=keycolor, $
;      xtickname=replicate(' ',10), ytickname=replicate(' ',10)
;    im_plothist, agesparent.z, bin=zbin, weight=agesparent.final_weight, $
;      /fill, fcolor=djs_icolor('orange'), /fline, /overplot, $
;      forientation=45, fspacing=fspacing, color=djs_icolor('orange')
;    im_plothist, agesmz.z, bin=zbin, weight=agesmz.final_weight, $
;      /overplot, thick=8.0, color=keycolor
;    plot, [0], [0], /nodata, /noerase, xtitle=xtitle, ytitle=htitle, $
;      xsty=1, ysty=1, xrange=xrange, yrange=hrange, position=hpos, $
;      color=keycolor, charsize=1.2

; redshift vs stellar mass
    plot, [0], [0], /nodata, /noerase, ysty=1, xsty=9, $
      xrange=xrange, yrange=yrange2, xtitle=xtitle, ytitle=ytitle2, $
      position=pos[*,1], yminor=5, color=keycolor
    ucsd_ages_oplot, agesmz.z, agesmz.k_mass, keynote=keynote

    im_plotconfig, /psclose
    if keyword_set(keynote) then begin
       spawn, 'ps2pdf13 '+psfile+' '+repstr(psfile,'.ps','.pdf'), /sh
       rmfile, psfile
    endif

stop    
    
; ---------------------------------------------------------------------------
; redshift histogram for AGES
    psfile = talkpath+'ages_zhist.ps'
    im_plotconfig, 0, pos, psfile=psfile, xmargin=[1.3,0.4], $
      width=6.8, height=4.5, keynote=keynote

    if keyword_set(keynote) then keycolor = djs_icolor('white') else $
      keycolor = djs_icolor('default')

    zbin = 0.02
    im_plothist, agesparent.z, bin=zbin, xbin, ybin, $
      weight=agesparent.final_weight, /noplot

    xtitle = 'Redshift'
    ytitle = 'Number of Galaxies'
    xrange = [0.0,0.8]
    yrange = [0,max(ybin)*1.1]

    if keyword_set(keynote) then $
      histcolor = fsc_color('wheat',101) else $
      histcolor = djs_icolor('default')
    
    plot, [0], [0], /nodata, xsty=5, ysty=5, xrange=xrange, yrange=yrange, $
      position=pos, xtickname=replicate(' ',10), ytickname=replicate(' ',10), $
      color=keycolor
    im_plothist, agesparent.z, bin=zbin, weight=agesparent.final_weight, $
      /overplot, /fill, fcolor=djs_icolor('orange'), /fline, $
      forientation=45, fspacing=fspacing, color=djs_icolor('orange')
;   im_plothist, sample.z, bin=zbin, weight=sample.final_weight, $
;     /overplot, /fill, fcolor=djs_icolor('orange'), /fline, $
;     forientation=135, fspacing=fspacing, color=djs_icolor('orange')
    im_plothist, agesmz.z, bin=zbin, weight=agesmz.final_weight, $
      /overplot, thick=8.0, color=histcolor

    plot, [0], [0], /nodata, /noerase, xtitle=xtitle, ytitle=ytitle, $
      xsty=1, ysty=1, xrange=xrange, yrange=yrange, position=pos, $
      color=keycolor
    legend, textoidl(['AGES I_{Vega}<20','Abundance Sample']), $
      /right, /top, box=0, textcolor=[djs_icolor('orange'),histcolor], $
      charsize=1.5, charthick=3.5, spacing=2.3, margin=0
    
    im_plotconfig, /psclose
    if keyword_set(keynote) then begin
       spawn, 'ps2pdf13 '+psfile+' '+repstr(psfile,'.ps','.pdf'), /sh
       rmfile, psfile
    endif

stop    
    
; ------------------------------------------------------------
; luminosity evolution for blue galaxies (literature compilation)
    psfile = talkpath+'redshift_vs_mg_lit.ps'
    im_plotconfig, 0, pos, psfile=psfile, charsize=2.0, $
      xmargin=[1.6,0.4], width=6.5, height=6.5, keynote=keynote

    if keyword_set(keynote) then keycolor = djs_icolor('white') else $
      keycolor = djs_icolor('default')

    h100 = 0.7
    B2g01 = +0.0759 ; ^{0.1}g = B+0.0759+0.0620*[(B-V)-0.5870] [AB, Blanton & Roweis]
    B2r01 = -0.6429 ; ^{0.1}r = B-0.6429-1.0845*[(B-V)-0.5870] [AB, Blanton & Roweis]
    Bvega2ab = (k_vega2ab(filterlist='bessell_B.par',/kurucz,/silent))[0]

; luminosity function colors from the literature
    if keyword_set(keynote) then begin
       w06color   = 'orange'       & w06sym   = 16 & w06psize = 2.0
       f06color   = 'dodger blue'  & f06sym   = 14 & f06psize = 3.0
       e07color   = 'tomato'       & e07sym   = 15 & e07psize = 2.0
       b06color   = 'dodger blue'  & b06sym   = 34 & b06psize = 3.0
    endif else begin
       w06color   = 'orange'       & w06sym   = 16 & w06psize = 2.0
       f06color   = 'navy'         & f06sym   = 14 & f06psize = 3.0
       e07color   = 'firebrick'    & e07sym   = 15 & e07psize = 2.0
       b06color   = 'navy'         & b06sym   = 34 & b06psize = 3.0
    endelse

    xtitle = 'Redshift'
    ytitle = textoidl('M_{0.1g}^{*} for Blue Galaxies') ; hack!
    xrange = [-0.02,0.85]
    yrange = [-20.5,-22.3]

; overplot the line for 1.5 mag/z of luminosity evolution    
    plot, [0], [0], /nodata, xtitle=xtitle, ytitle=ytitle, $
      xsty=1, ysty=1, xrange=xrange, yrange=yrange, position=pos, $
      color=keycolor
    djs_oplot, zaxis1, poly(zaxis1-0.1,[-20.95,-1.5]), line=0, color=keycolor
    im_legend, 'Q=1.5 mag z^{-1}', /right, /bottom, box=0, $
      charsize=1.7, line=0, color=keycolor, textcolor=keycolor, $
      pspacing=1.9, thick=6.0
    
; -------------------------
; Eisenstein et al. 2009 [AGES]; Omega_0=0.3, Omega_lamba=0.7, h=1.0;
; AB; alpha fixed at -1.10 for all blue galaxies; evolving color cut
; based on the (u0.1-r0.1) color: A=u0.1-r0.1+0.08(M_r0.1+20);  

; updated values from Cool+10!
    z_e07    = [0.1,0.2,0.3,0.4,0.5,0.65]
    zerr_e07 = [0.05,0.05,0.05,0.05,0.05,0.1]

    mstar_e07 = [-20.03,-20.42,-20.51,-20.72,-20.93,-20.90] + 5.0*alog10(h100) ; h=1-->h=0.7
    mstarerr_e07 = [0.12,0.07,0.04,0.05,0.11,0.08]
;   mstar_e07 = [-19.95,-20.39,-20.43,-20.57,-20.86,-20.83] + 5.0*alog10(h100) ; h=1-->h=0.7
;   mstarerr_e07 = [0.12,0.07,0.04,0.05,0.11,0.08]

; -------------------------
; Willmer et al. 2006 [DEEP2]; "minimal" model adopted, as
; recommended; Omega_0=0.3, Omega_lamba=0.7, h=0.7; Vega; alpha=-1.30
; fixed for all blue galaxies, based on a fit to all galaxies at
; 0.2<z<0.6 in COMBO-17 (Faber et al. 2005); color cut: (U-B) =
; -0.032(M_B+21.52)+0.454-0.25; drop the last bin
    z_w06 = [0.3,0.5,0.7,0.9,1.1];,1.3]
    zerr_w06 = [0.1,0.1,0.1,0.1,0.1];,0.1]

    mstar_w06 = [-20.36,-20.72,-21.15,-21.21,-21.38] + Bvega2ab + B2r01 ; Vega-->AB; B-->r0.1 ; ,-21.86]
    mstarerr_up_w06 = [0.13,0.05,0.07,0.01,0.04];,0.07]
    mstarerr_lo_w06 = [0.11,0.07,0.07,0.03,0.05];,0.08]
    mstarerr_w06 = mstar_w06*0.0
    for i = 0, n_elements(mstar_w06)-1L do mstarerr_w06[i] = mean([mstarerr_up_w06[i],mstarerr_lo_w06[i]])
    
; -------------------------
; Faber et al. 2007 [COMBO17]; Omega_0=0.3, Omega_lamba=0.7, h=0.7;
; Vega; alpha=-1.3 fixed for all blue galaxies;
    z_f06 = [0.3,0.5,0.7,0.9,1.1]
    zerr_f06 = [0.1,0.1,0.1,0.1,0.1]

    mstar_f06 = [-20.74,-21.10,-21.30,-21.10,-21.25] + Bvega2ab + B2r01 ; Vega-->AB; B-->r0.1
    mstarerr_f06 = [0.20,0.15,0.16,0.17,0.18]

    oploterror, z_f06, mstar_f06, zerr_f06, mstarerr_f06, symsize=3.8, $
      psym=symcat(f06sym), color=fsc_color(f06color,1E8), errcolor=fsc_color(f06color,1E8)
    
    oploterror, z_w06, mstar_w06, zerr_w06, mstarerr_up_w06, /hi, symsize=3.5, $
      psym=symcat(w06sym), color=fsc_color(w06color,1E8), errcolor=fsc_color(w06color,1E8)
    oploterror, z_w06, mstar_w06, zerr_w06, mstarerr_lo_w06, /lo, symsize=3.5, $
      psym=symcat(w06sym), color=fsc_color(w06color,1E8), errcolor=fsc_color(w06color,1E8)

    oploterror, z_e07, mstar_e07, zerr_e07, mstarerr_e07, symsize=3.5, $
      psym=symcat(e07sym), color=fsc_color(e07color,1E8), errcolor=fsc_color(e07color,1E8)

; legend
;   label = ['AGES+SDSS','COMBO-17','DEEP2']
    label = ['Willmer+06 (DEEP2)','Faber+07 (COMBO17)','Cool+10 (AGES)']
    im_legend, label, /left, /top, box=0, charsize=1.5, $
      psym=[w06sym,f06sym,e07sym], symsize=[2.4,2.8,2.0], $
      spacing=2.0, color=[w06color,f06color,e07color], $
      textcolor=[w06color,f06color,e07color]

    im_plotconfig, /psclose
    if keyword_set(keynote) then begin
       spawn, 'ps2pdf13 '+psfile+' '+repstr(psfile,'.ps','.pdf'), /sh
       rmfile, psfile
    endif

stop    
    
return
end
