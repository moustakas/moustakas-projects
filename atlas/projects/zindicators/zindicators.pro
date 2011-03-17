;+
; NAME:
;       ZINDICATORS
;
; CALLING SEQUENCE:
;
; PURPOSE:
;       Investigate the electron-temperature dependence of various
;       line ratios in HII regions.
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2004-2006, U of A
;       jm08may13nyu - revisited
;-

pro zindicators, hii, postscript=postscript, pdf=pdf, $
  encapsulated=encapsulated, _extra=extra

; path names    

    abundpath = atlas_path(/projects)+'abundances/'
    paperpath = atlas_path(/papers)+'abundances/'
    pspath = paperpath+'FIG_ABUNDANCES/'
    
; read the data and the models    
    
    if (n_elements(hii) eq 0L) then hii = read_hii_regions(/nosdss,_extra=extra)

    sbgrids = read_kewley_grids(model=3,Z=Z,U=U)  ; model=1 (SB, Pegase, continuous)
    hiigrids = read_kewley_grids(model=8,Z=Z,U=U) ; model=6 (HII, Pegase, instantaneous)

; initialize the plotting variables

    postthick1 = 2.0
    postthick2 = 2.0
    postthick3 = 2.0
    postthick4 = 2.0
    postthick5 = 2.0
    postthick6 = 1.0
    postthick7 = 2.0

    if keyword_set(pdf) then begin

       postscript = 0L
       encapsulated = 0L
       pspath = paperpath+'/keynote/'

       textcolor1 = 'white'
       axis_color = 'black'

       postthick1 = 6.0 
       postthick2 = 6.0 
       postthick3 = 8.0 
       postthick4 = 4.0 
       postthick5 = 8.0 
       postthick6 = 1.0 
       postthick7 = 10.0 

    endif
    
    if keyword_set(postscript) then begin
    
       textcolor1 = 'black'
       axis_color = 'black'
 
       postthick1 = 4.0 
       postthick2 = 3.0 
       postthick3 = 8.0 
       postthick4 = 4.0 
       postthick5 = 6.0 
       postthick6 = 1.0 
       postthick7 = 10.0 

    endif

    if keyword_set(encapsulated) then armexten = '.eps' else armexten = '.ps'

    if keyword_set(pdf) then begin
       hiipsize1 = 0.6 & hiisym1 = 108 & hiicolor1 = textcolor1
    endif else begin
       hiipsize1 = 0.6 & hiisym1 = 108 & hiicolor1 = 'dark grey'
    endelse
    
    @'xyrange_zindicators'

; additional constants    
    
;   zsun = 0.02 & log12ohsun = 8.93
    zsun = 0.0122 & log12ohsun = 8.66
    log12ohgrid = alog10(Z*0.02/zsun)+log12ohsun
    dlog12oh = 0.005

    log12ohplot = findgen((max(log12ohgrid)-min(log12ohgrid))/dlog12oh)*dlog12oh+min(log12ohgrid)
    nlog12oh = n_elements(log12ohplot)

    logu1 = 0L & logu2 = 5L
    
; ###########################################################################    
; Paper Plots
; ###########################################################################    
    
; ------------------------------------------------------------
; 3-panel [N II]/Ha, ([O III]/Hb)/([N II]/Ha), and R23 vs 12+log(O/H) [Te]
; ------------------------------------------------------------

    psname = 'hii_ratios_vs_12oh_te_3panel'
    if keyword_set(postscript) or keyword_set(pdf) then begin
       psfile = pspath+psname+armexten
       if keyword_set(postscript) then splog, 'Writing '+psfile
    endif else delvarx, psfile

    xmargin = [1.1,0.4] & ymargin = [0.4,1.1]
    xspace = 0.0 & yspace = 0.0
    width = 3.0*[1,1,1] & height = 3.0
    xpage = total(height)+total(ymargin)+total(yspace)
    ypage = total(width)+total(xmargin)+total(xspace)

    xspace1 = xspace & yspace1 = yspace
    xmargin1 = xmargin & ymargin1 = ymargin
    width1 = width & height1 = height
    xpage1 = xpage & ypage1 = ypage
    
    arm_plotconfig, landscape=1, nx=3, ny=1, xmargin=xmargin1, $
      ymargin=ymargin1, xspace=xspace1, yspace=yspace1, width=width1, $
      height=height1, coord=pos, xpage=xpage1, ypage=ypage1, $
      psfile=psfile, /writeover, bw=0L;, /show
    cleanplot, /silent

    im_symbols, hiisym1, psize=hiipsize1, thick=postthick1, fill=1, color=fsc_color(hiicolor1,150)

    yrange = [6.85,9.35]
    ytitle = '12 + log (O/H) [T_{e}]'
    
; (a) [N II]/Ha; van Zee et al. (1998), Denicolo et al. (2002),
; Storchi-Bergmann, Calzetti, & Kinney (1994), Raimann et al. (2000),
; Pettini & Pagel (2004)
    
    indx = where((hii.zstrong_niiha gt -900.0) and (hii.zt_log12oh_te gt -900.0),nindx)
    
    x = hii[indx].zstrong_niiha
    xerr = hii[indx].zstrong_niiha_err

    y = hii[indx].zt_log12oh_te
    yerr = hii[indx].zt_log12oh_te_err

    xrange = [-3,0.4]
;   xtitle = 'log (N2)'
    xtitle = 'log ([N II]/H\alpha)'
    
    plot, [0], [0], /nodata, xsty=1, ysty=1, xrange=xrange, yrange=yrange, $
      xthick=postthick1, ythick=postthick1, charthick=postthick2, charsize=charsize_3, $
      xtitle=textoidl(xtitle), ytitle=textoidl(ytitle), color=fsc_color(textcolor1,100), $
      position=pos[*,0], xtickinterval=1.0
    oplot, x, y, psym=8, thick=postthick1, color=fsc_color(hiicolor1,150)
;   legend, '(a)', /left, /top, box=0, charsize=charsize_4, charthick=postthick

; overplot various calibrations

    niiha_denicolo = findgen((0.0-(-2.6))/0.01+1)*0.01+(-2.6)
    ohaxis_denicolo = 0.73*niiha_denicolo + 9.12
    niiha_melbourne = findgen((-0.45-min(x))/0.01+1)*0.01+min(x)
    ohaxis_melbourne = 9.26 + 1.23*niiha_melbourne + 0.204*niiha_melbourne^2
;   ohaxis_vanzee = findgen((9.1-7.0)/0.01+1)*0.01+7.0
;   niiha_vanzee = (ohaxis_vanzee-9.36) / 1.02
    ohaxis_pettini = findgen((8.729-7.475)/0.01+1)*0.01+7.475 ; only valid in this range!
    niiha_pettini = (ohaxis_pettini-8.90) / 0.57

    oplot, niiha_denicolo, ohaxis_denicolo, line=5, thick=postthick3, color=fsc_color('orange',194)
    oplot, niiha_pettini, ohaxis_pettini, line=0, thick=postthick3, color=fsc_color('orchid',189)
;   djs_oplot, niiha_melbourne, ohaxis_melbourne, line=5, thick=postthick, color='orange'

; fit a linear bisector 

    niihamin = min(x) & niihamax = max(x) ; -0.5
    
    range = where((x ge niihamin) and (x le niihamax))
    fitx = x[range] & fity = y[range] & fitxerr = xerr[range] & fityerr = yerr[range]

    sixlin, fitx, fity, a, siga, b, sigb
    coeff = [a[2],b[2]] & coeff_err = [siga[2],sigb[2]]

    xaxis = findgen((max(fitx)-min(fitx))/0.001+1)*0.001+min(fitx)
    yfit = poly(xaxis,coeff)
;   djs_oplot, xaxis, yfit, line=0, thick=postthick2;, color='green'

    splog, '[N II]/Ha vs 12+log(O/H) [Te] coefficients: '
    niceprint, coeff, coeff_err, minmax(fitx)

; ([O III]/Hb)/([N II]/Ha) vs 12+log (O/H) [Te] - Pettini & Pagel (2004)

    indx = where((hii.oiii_5007_h_beta_nii_6584_h_alpha gt -900.0) and (hii.zt_log12oh_te gt -900.0),nindx)
    
    x = hii[indx].oiii_5007_h_beta_nii_6584_h_alpha
    xerr = hii[indx].oiii_5007_h_beta_nii_6584_h_alpha_err

    y = hii[indx].zt_log12oh_te
    yerr = hii[indx].zt_log12oh_te_err

    xrange = [-1.2,3.8]
;   xtitle = 'log (O3N2)'
    xtitle = 'log {([O III]/H\beta)/([N II]/H\alpha)}'
    
    plot, [0], [0], /nodata, xsty=1, ysty=1, xrange=xrange, yrange=yrange, $
      xthick=postthick1, ythick=postthick1, charthick=postthick2, charsize=charsize_3, $
      xtitle=textoidl(xtitle), ytitle=textoidl(''), color=fsc_color(textcolor1,100), $
      position=pos[*,1], /noerase, ytickname=replicate(' ',10)
    oplot, x, y, psym=8, thick=postthick1, color=fsc_color(hiicolor1,150)
;   legend, '(b)', /left, /top, box=0, charsize=charsize_4, charthick=postthick

; overplot the Pettini & Pagel (2004) relationship; only applicable
; from -1.0 to 1.9 in oiiinii

    o3n2_pettini = findgen((1.9-(-1.0))/0.01+1)*0.01+(-1.0)
    oh_pettini = -0.32*o3n2_pettini + 8.73 
    oplot, o3n2_pettini, oh_pettini, line=0, thick=postthick3, color=fsc_color('orchid',189)

; fit a linear function/bisector

    range = where((x ge min(x)) and (x le 1.9))
    fitx = x[range]
    fity = y[range]
    fitxerr = xerr[range]
    fityerr = yerr[range]

    sixlin, fitx, fity, a, siga, b, sigb
    coeff = [a[2],b[2]] & coeff_err = [siga[2],sigb[2]]

; note I changed min(fitx) to -1.0!

    xaxis = findgen(((1.9)-(-1.0))/0.001+1)*0.001+(-1.0)
;   xaxis = findgen((max(fitx)-min(fitx))/0.001+1)*0.001+min(fitx)
    yfit = poly(xaxis,coeff)
;   djs_oplot, xaxis, yfit, line=0, thick=postthick2;, color='green'

    splog, 'oiiinii vs 12+log(O/H) [Te] coefficients and range: '
    niceprint, coeff, coeff_err, minmax(fitx)

; R23 vs 12+log (O/H) [Te]

    indx = where((hii.zstrong_r23 gt -900.0) and (hii.zt_log12oh_te gt -900.0),nindx)
     
    x = alog10(hii[indx].zstrong_r23)
    xerr = hii[indx].zstrong_r23_err/hii[indx].zstrong_r23/alog(10.0)

    y = hii[indx].zt_log12oh_te
    yerr = hii[indx].zt_log12oh_te_err

    xrange = [-0.49,1.2]
;   xtitle = 'log (R_{23})'
    xtitle = 'log {([O II]+[O III])/H\beta}'
    
    plot, [0], [0], /nodata, xsty=1, ysty=1, xrange=xrange, yrange=yrange, $
      xthick=postthick1, ythick=postthick1, charthick=postthick2, charsize=charsize_3, $
      xtitle=textoidl(xtitle), ytitle=textoidl(''), color=fsc_color(textcolor1,100), $
      position=pos[*,2], /noerase, ytickname=replicate(' ',10)
    oplot, x, y, psym=8, thick=postthick1, color=fsc_color(hiicolor1,150)
;   legend, '(c)', /left, /top, box=0, charsize=charsize_4, charthick=postthick

; overplot the Kewley and McGaugh theoretical calibrations
 
    logo32 = 0.0
    logq = alog10(4E7)

    r23max = 0.97 & r23min = -0.2 & dr23 = 0.01
    logr23 = findgen((r23max-r23min)/dr23+1)*dr23+r23min
 
    M91_upper = 12.0 - 2.939 - 0.2*logr23 - 0.237*logr23^2 - 0.305*logr23^3 - $
      0.0283*logr23^4 - logo32*(0.0047 - 0.0221*logr23 - 0.102*logr23^2 - $
      0.0817*logr23^3 - 0.00717*logr23^4)
    m91_lower = (12.0 - 4.944 + 0.767*logr23 + 0.602*logr23^2) - logo32*(0.29 + 0.332*logr23 - 0.331*logr23^2)

    KK04_upper = 9.72D - 0.777*logr23 - 0.951*logr23^2 - 0.072*logr23^3 - 0.811*logr23^4 - $
      logq*(0.0737 - 0.0713*logr23 - 0.141*logr23^2 + 0.0373*logr23^3 - 0.058*logr23^4)
    kk04_lower = 9.40D + 4.65D*logr23 - 3.17D*logr23^2 - logq*(0.272D + 0.547D*logr23 - 0.513D*logr23^2)

    oplot, logr23, M91_upper, line=0, thick=postthick3, color=fsc_color('forest green',190)
    oplot, logr23, M91_lower, line=0, thick=postthick3, color=fsc_color('forest green',190)
    oplot, logr23, KK04_upper, line=2, thick=postthick3, color=fsc_color('dodger blue',191)
    oplot, logr23, KK04_lower, line=2, thick=postthick3, color=fsc_color('dodger blue',191)

; overplot the Pilyugin & Thuan empirical calibration

    p = 0.72
    pt05_lower = (10^logr23 + 106.4 + 106.8*P - 3.40*P^2) / (17.72 + 6.60*P + 6.95*P^2 - 0.302*10^logr23)
    pt05_upper = (10^logr23 + 726.1 + 842.2*P + 337.5*P^2) / (85.96 + 82.76*P + 43.98*P^2 + 1.793*10^logr23)

;   djs_oplot, logr23, pt05_upper, line=3, thick=postthick3, color=fsc_color('firebrick',192)
;   djs_oplot, logr23, pt05_lower, line=3, thick=postthick3, color=fsc_color('firebrick',192)
    
    im_openclose, postscript=(keyword_set(postscript) or keyword_set(pdf)), /close    
    if keyword_set(pdf) then begin
       splog, 'Writing '+pspath+psname+'.pdf'
       spawn, 'ps2pdf13 '+pspath+psname+'.ps '+pspath+psname+'.pdf', /sh
       rmfile, pspath+psname+'.ps' ; delete the junk ps file
    endif

stop    

; ------------------------------------------------------------
; R23 vs 12+log(O/H) [Te]
; ------------------------------------------------------------

    psname = 'hii_12oh_te_vs_r23'
    xpage = 8.5 & ypage = 8.0
    im_openclose, pspath+psname, postscript=(keyword_set(postscript) or $
      keyword_set(pdf)), xsize=xpage, ysize=ypage, encapsulated=encapsulated, $
      silent=keyword_set(pdf)

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.6, height=6.6, $
      xmargin=[1.5,0.4], ymargin=[0.3,1.1], xpage=8.5, ypage=8.0, $
      position=pos, /normal

    indx = where((hii.zstrong_r23 gt -900.0) and (hii.zt_log12oh_te gt -900.0),nindx)
     
    x = alog10(hii[indx].zstrong_r23)
    xerr = hii[indx].zstrong_r23_err/hii[indx].zstrong_r23/alog(10.0)

    y = hii[indx].zt_log12oh_te
    yerr = hii[indx].zt_log12oh_te_err

    xrange = r23range3
    yrange = ohrange2
    
    xtitle = 'log {([O II] + [O III]) / H\beta}'
;   xtitle = 'log (R_{23}) = ([O II] + [O III]) / H\beta'
    ytitle = '12 + log (O/H) [T_{e}]'

    plot, [0], [0], /nodata, xsty=1, ysty=1, xrange=xrange, yrange=yrange, $
      xthick=postthick1, ythick=postthick1, charthick=postthick2, charsize=charsize_8, $
      xtitle=textoidl(xtitle), ytitle=textoidl(ytitle), color=fsc_color(textcolor1,100), $
      position=pos[*,0]
    im_symbols, hiisym1, psize=hiipsize1, thick=postthick1, fill=1, color=fsc_color(hiicolor1,150)
    oplot, x, y, psym=8, thick=postthick1, color=fsc_color(hiicolor1,150)
;   oploterror, x, y, xerr, yerr, psym=8, thick=postthick1, color=fsc_color(hiicolor1,150), $
;     errthick=postthick1, errcolor=fsc_color(hiicolor1,150)
    
; overplot the Kewley and McGaugh theoretical calibrations
 
    logo32 = 0.0
    logq = alog10(4E7)

    r23max = 0.97 & r23min = -0.2 & dr23 = 0.01
    logr23 = findgen((r23max-r23min)/dr23+1)*dr23+r23min
 
    M91_upper = 12.0 - 2.939 - 0.2*logr23 - 0.237*logr23^2 - 0.305*logr23^3 - $
      0.0283*logr23^4 - logo32*(0.0047 - 0.0221*logr23 - 0.102*logr23^2 - $
      0.0817*logr23^3 - 0.00717*logr23^4)
    m91_lower = (12.0 - 4.944 + 0.767*logr23 + 0.602*logr23^2) - logo32*(0.29 + 0.332*logr23 - 0.331*logr23^2)
    m91_keep = where(m91_upper gt m91_lower)
    m91_upper = m91_upper[m91_keep] & m91_lower = m91_lower[m91_keep]

    KK04_upper = 9.72D - 0.777*logr23 - 0.951*logr23^2 - 0.072*logr23^3 - 0.811*logr23^4 - $
      logq*(0.0737 - 0.0713*logr23 - 0.141*logr23^2 + 0.0373*logr23^3 - 0.058*logr23^4)
    kk04_lower = 9.40D + 4.65D*logr23 - 3.17D*logr23^2 - logq*(0.272D + 0.547D*logr23 - 0.513D*logr23^2)
    kk04_keep = where(kk04_upper gt kk04_lower)
    kk04_upper = kk04_upper[kk04_keep] & kk04_lower = kk04_lower[kk04_keep]
    
    djs_oplot, logr23, M91_upper, line=0, thick=postthick3, color=fsc_color('forest green',190)
    djs_oplot, logr23, M91_lower, line=0, thick=postthick3, color=fsc_color('forest green',190)
    djs_oplot, logr23, KK04_upper, line=5, thick=postthick3, color=fsc_color('dodger blue',191)
    djs_oplot, logr23, KK04_lower, line=5, thick=postthick3, color=fsc_color('dodger blue',191)

; overplot the Pilyugin & Thuan empirical calibration

    p = 0.72
    pt05_lower = (10^logr23 + 106.4 + 106.8*P - 3.40*P^2) / (17.72 + 6.60*P + 6.95*P^2 - 0.302*10^logr23)
    pt05_upper = (10^logr23 + 726.1 + 842.2*P + 337.5*P^2) / (85.96 + 82.76*P + 43.98*P^2 + 1.793*10^logr23)

    djs_oplot, logr23, pt05_upper, line=3, thick=postthick3, color=fsc_color('firebrick',192)
    djs_oplot, logr23, pt05_lower, line=3, thick=postthick3, color=fsc_color('firebrick',192)

    im_openclose, postscript=(keyword_set(postscript) or keyword_set(pdf)), /close    
    if keyword_set(pdf) then begin
       splog, 'Writing '+pspath+psname+'.pdf'
       spawn, 'ps2pdf13 '+pspath+psname+'.ps '+pspath+psname+'.pdf', /sh
       rmfile, pspath+psname+'.ps' ; delete the junk ps file
    endif

; ------------------------------------------------------------
; 4-panel 12+log(O/H) [Te] vs [OIII]-based line-ratios
; ------------------------------------------------------------

    psname = 'hii_12oh_te_vs_oiiiratios_4panel'
    im_openclose, pspath+psname, postscript=postscript, encapsulated=encapsulated, xsize=8.5, ysize=7.4

    pagemaker, nx=2, ny=2, xspace=0.0, yspace=0.0, width=3.05, height=3.05, $
      xmargin=[1.3,1.1], ymargin=[0.2,1.1], xpage=8.5, ypage=7.4, $
      position=pos, /normal
    
    xtitle = '12 + log (O/H) [T_{e}]'
    xrange = ohrange

; --------------------------------------------------    
; Panel 1 - [O III]/Hb
; --------------------------------------------------    

    indx = where((hii.zstrong_oiiihb gt -900.0) and (hii.zt_log12oh_te gt -900.0),nindx)
    
    x = hii[indx].zt_log12oh_te
    xerr = hii[indx].zt_log12oh_te_err

    y = hii[indx].zstrong_oiiihb
    yerr = hii[indx].zstrong_oiiihb_err

    yrange = oiiihbrange
    ytitle = 'log ([O III]/H\beta)'

    zindicators_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, talk=dotalk, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=3, ystyle=3, /right, /top, position=pos[*,0], xtickname=replicate(' ',10), $
      charsize=charsize_5

; overlay the kewley grids

    hiimodel1 = interpol(hiigrids[*,logu1].oiii_5007_h_beta,log12ohgrid,log12ohplot,/quad)
    hiimodel2 = interpol(hiigrids[*,logu2].oiii_5007_h_beta,log12ohgrid,log12ohplot,/quad)
    sbmodel1 = interpol(sbgrids[*,logu1].oiii_5007_h_beta,log12ohgrid,log12ohplot,/quad)
    sbmodel2 = interpol(sbgrids[*,logu2].oiii_5007_h_beta,log12ohgrid,log12ohplot,/quad)    

    djs_oplot, log12ohplot, hiimodel1, line=0, thick=postthick, color=djs_icolor(talkcolor)
    djs_oplot, log12ohplot, hiimodel2, line=0, thick=postthick, color=djs_icolor(talkcolor)
;   djs_oplot, log12ohplot, sbmodel1, line=0, thick=postthick, color=djs_icolor(talkcolor)
;   djs_oplot, log12ohplot, sbmodel2, line=0, thick=postthick, color=djs_icolor(talkcolor)
    
    xyouts, log12ohplot[logu1], hiimodel1[logu1], string(U[logu1],format='(F5.2)'), /data, charsize=1.1, $
      charthick=postthick, align=1.0, color=djs_icolor(talkcolor)
    xyouts, log12ohplot[logu1], hiimodel2[logu1], string(U[logu2],format='(F5.2)'), /data, charsize=1.1, $
      charthick=postthick, align=1.0, color=djs_icolor(talkcolor)

;   legend, ['HII Region','Starburst'], line=[2,0], /right, /bottom, box=0, $
;     charsize=1.5, charthick=postthick, thick=postthick
       
    legend, '(a)', /left, /top, box=0, charsize=charsize_5, charthick=postthick, $
      textcolor=djs_icolor(talkcolor)

; --------------------------------------------------    
; Panel 2 - [O II]/Hb
; --------------------------------------------------    

    indx = where((hii.oii_h_beta gt -900.0) and (hii.zt_log12oh_te gt -900.0),nindx)
    
    x = hii[indx].zt_log12oh_te
    xerr = hii[indx].zt_log12oh_te_err

    y = hii[indx].oii_h_beta
    yerr = hii[indx].oii_h_beta_err

    yrange = oiihbrange
    ytitle = 'log ([O II]/H\beta)'
    
    zindicators_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, talk=dotalk, $
      xtitle='', ytitle='', ysty=11, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=3, /right, /top, position=pos[*,1], /noerase, ytickname=replicate(' ',10), $
      xtickname=replicate(' ',10), charsize=charsize_5, hiipsize=psize, /errorleft
    axis, /yaxis, ythick=postthick, ysty=1, yrange=!y.crange, ytitle=textoidl(ytitle), $
      charthick=postthick, charsize=charsize_5, color=djs_icolor(talkcolor)

; overlay the kewley grids

    hiimodel1 = interpol(hiigrids[*,logu1].oii_h_beta,log12ohgrid,log12ohplot,/quad)
    hiimodel2 = interpol(hiigrids[*,logu2].oii_h_beta,log12ohgrid,log12ohplot,/quad)
    sbmodel1 = interpol(sbgrids[*,logu1].oii_h_beta,log12ohgrid,log12ohplot,/quad)
    sbmodel2 = interpol(sbgrids[*,logu2].oii_h_beta,log12ohgrid,log12ohplot,/quad)
    
    djs_oplot, log12ohplot, hiimodel1, line=0, thick=postthick, color=djs_icolor(talkcolor)
    djs_oplot, log12ohplot, hiimodel2, line=0, thick=postthick, color=djs_icolor(talkcolor)
;   djs_oplot, log12ohplot, sbmodel1, line=0, thick=postthick
;   djs_oplot, log12ohplot, sbmodel2, line=0, thick=postthick
    
    xyouts, log12ohplot[nlog12oh-1], hiimodel1[nlog12oh-1], string(U[logu1],format='(F5.2)'), /data, charsize=1.1, $
      charthick=postthick, align=0.0, color=djs_icolor(talkcolor)
    xyouts, log12ohplot[nlog12oh-1], hiimodel2[nlog12oh-1], string(U[logu2],format='(F5.2)'), /data, charsize=1.1, $
      charthick=postthick, align=0.0, color=djs_icolor(talkcolor)
;   xyouts, log12ohplot[logu1], hiimodel1[logu1], string(U[logu1],format='(F5.2)'), /data, charsize=1.1, charthick=postthick, align=1.0
;   xyouts, log12ohplot[logu1], hiimodel2[logu1], string(U[logu2],format='(F5.2)'), /data, charsize=1.1, charthick=postthick, align=1.0

;   legend, ['HII Region','Starburst'], line=[2,0], /left, /bottom, box=0, $
;     charsize=1.5, charthick=postthick, thick=postthick

    legend, '(b)', /left, /top, box=0, charsize=charsize_5, charthick=postthick, $
      textcolor=djs_icolor(talkcolor)

; --------------------------------------------------    
; Panel 3 - O32
; --------------------------------------------------    
 
    indx = where((hii.zstrong_o32 gt -900.0) and (hii.zt_log12oh_te gt -900.0),nindx)
    
    x = hii[indx].zt_log12oh_te
    xerr = hii[indx].zt_log12oh_te_err

    y = hii[indx].zstrong_o32
    yerr = hii[indx].zstrong_o32_err

    yrange = o32range ; [-3,0]
    ytitle = 'log (O_{32})'
;   ytitle = 'log ([O III]/[O II])'
    
    zindicators_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, talk=dotalk, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=3, ystyle=3, /right, /top, position=pos[*,2], /noerase, $
      charsize=charsize_5, hiipsize=psize
 
; overlay some representative kewley grids

    hiimodel1 = interpol(hiigrids[*,logu1].o32,log12ohgrid,log12ohplot,/quad)
    hiimodel2 = interpol(hiigrids[*,logu2].o32,log12ohgrid,log12ohplot,/quad)
    sbmodel1 = interpol(sbgrids[*,logu1].o32,log12ohgrid,log12ohplot,/quad)
    sbmodel2 = interpol(sbgrids[*,logu2].o32,log12ohgrid,log12ohplot,/quad)
    
    djs_oplot, log12ohplot, hiimodel1, line=0, thick=postthick, color=djs_icolor(talkcolor)
    djs_oplot, log12ohplot, hiimodel2, line=0, thick=postthick, color=djs_icolor(talkcolor)
;   djs_oplot, log12ohplot, sbmodel1, line=0, thick=postthick, color=djs_icolor(talkcolor)
;   djs_oplot, log12ohplot, sbmodel2, line=0, thick=postthick, color=djs_icolor(talkcolor)
    
    xyouts, log12ohplot[nlog12oh-1], hiimodel1[nlog12oh-1], string(U[logu1],format='(F5.2)'), /data, charsize=1.1, $
      charthick=postthick, align=0.0, color=djs_icolor(talkcolor)
    xyouts, log12ohplot[nlog12oh-1], hiimodel2[nlog12oh-1], string(U[logu2],format='(F5.2)'), /data, charsize=1.1, $
      charthick=postthick, align=0.0, color=djs_icolor(talkcolor)
;   xyouts, log12ohplot[logu1], hiimodel1[logu1], string(U[logu1],format='(F5.2)'), /data, charsize=1.1, charthick=postthick, align=1.0
;   xyouts, log12ohplot[logu1], hiimodel2[logu1], string(U[logu2],format='(F5.2)'), /data, charsize=1.1, charthick=postthick, align=1.0
   
;   legend, ['HII Region','Starburst'], line=[2,0], /right, /bottom, box=0, $
;     charsize=1.5, charthick=postthick, thick=postthick

    legend, '(c)', /left, /top, box=0, charsize=charsize_5, charthick=postthick, $
      textcolor=djs_icolor(talkcolor)

; --------------------------------------------------    
; Panel 4 - R23
; --------------------------------------------------    

    indx = where((hii.zstrong_r23 gt -900.0) and (hii.zt_log12oh_te gt -900.0),nindx)
    
    x = hii[indx].zt_log12oh_te
    xerr = hii[indx].zt_log12oh_te_err

    y = alog10(hii[indx].zstrong_r23)
    yerr = hii[indx].zstrong_r23_err/hii[indx].zstrong_r23/alog(10.0)

    yrange = r23range3
    ytitle = 'log (R_{23})'
    
    zindicators_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, talk=dotalk, $
      xtitle=xtitle, ytitle='', ystyle=11, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=3, /right, /top, position=pos[*,3], /noerase, $
      ytickname=replicate(' ',10), charsize=charsize_5, hiipsize=psize, /errorleft
    axis, /yaxis, ythick=postthick, ysty=1, yrange=!y.crange, ytitle=textoidl(ytitle), $
      charthick=postthick, charsize=charsize_5, color=djs_icolor(talkcolor)

; overlay some representative kewley grids

    hiimodel1 = interpol(hiigrids[*,logu1].r23,log12ohgrid,log12ohplot,/quad)
    hiimodel2 = interpol(hiigrids[*,logu2].r23,log12ohgrid,log12ohplot,/quad)
    sbmodel1 = interpol(sbgrids[*,logu1].r23,log12ohgrid,log12ohplot,/quad)
    sbmodel2 = interpol(sbgrids[*,logu2].r23,log12ohgrid,log12ohplot,/quad)
    
    djs_oplot, log12ohplot, hiimodel1, line=0, thick=postthick, color=djs_icolor(talkcolor)
    djs_oplot, log12ohplot, hiimodel2, line=0, thick=postthick, color=djs_icolor(talkcolor)
;   djs_oplot, log12ohplot, sbmodel1, line=0, thick=postthick, color=djs_icolor(talkcolor)
;   djs_oplot, log12ohplot, sbmodel2, line=0, thick=postthick, color=djs_icolor(talkcolor)
    
    xyouts, log12ohplot[nlog12oh-1], hiimodel1[nlog12oh-1], string(U[logu1],format='(F5.2)'), /data, charsize=1.1, $
      charthick=postthick, align=0.0, color=djs_icolor(talkcolor)
    xyouts, log12ohplot[nlog12oh-1], hiimodel2[nlog12oh-1], string(U[logu2],format='(F5.2)'), /data, charsize=1.1, $
      charthick=postthick, align=0.0, color=djs_icolor(talkcolor)
;   xyouts, log12ohplot[logu1], hiimodel1[logu1], string(U[logu1],format='(F5.2)'), /data, charsize=1.1, charthick=postthick, align=1.0
;   xyouts, log12ohplot[logu1], hiimodel2[logu1], string(U[logu2],format='(F5.2)'), /data, charsize=1.1, charthick=postthick, align=1.0

;   legend, ['HII Region','Starburst'], line=[2,0], /right, /top, box=0, $
;     charsize=1.5, charthick=postthick, thick=postthick

    legend, '(d)', /left, /top, box=0, charsize=charsize_5, charthick=postthick, $
      textcolor=djs_icolor(talkcolor)

    im_openclose, postscript=(keyword_set(postscript) or keyword_set(pdf)), /close    
    if keyword_set(pdf) then begin
       splog, 'Writing '+pspath+psname+'.pdf'
       spawn, 'ps2pdf13 '+pspath+psname+'.ps '+pspath+psname+'.pdf', /sh
       rmfile, pspath+psname+'.ps' ; delete the junk ps file
    endif

; ------------------------------------------------------------
; 4-panel 12+log(O/H) [Te] vs [NII]-based line-ratios
; ------------------------------------------------------------

    psname = 'hii_12oh_te_vs_niiratios_4panel'
    im_openclose, pspath+psname, postscript=postscript, encapsulated=encapsulated, xsize=8.5, ysize=7.4

    pagemaker, nx=2, ny=2, xspace=0.0, yspace=0.0, width=3.05, height=3.05, $
      xmargin=[1.3,1.1], ymargin=[0.2,1.1], xpage=8.5, ypage=7.4, $
      position=pos, /normal
    
    xtitle = '12 + log (O/H) [T_{e}]'
    xrange = ohrange

; --------------------------------------------------    
; Panel 1 - [N II]/[O II]
; --------------------------------------------------    

    indx = where((hii.zstrong_niioii gt -900.0) and (hii.zt_log12oh_te gt -900.0),nindx)
    
    x = hii[indx].zt_log12oh_te
    xerr = hii[indx].zt_log12oh_te_err

    y = hii[indx].zstrong_niioii
    yerr = hii[indx].zstrong_niioii_err

    yrange = niioiirange
    ytitle = 'log ([N II]/[O II])'
;   ytitle = 'log ([N II] \lambda6584 / [O II] \lambda3727)'

    zindicators_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, talk=dotalk, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=3, ystyle=3, /right, /top, position=pos[*,0], xtickname=replicate(' ',10), $
      charsize=charsize_5

; overlay the kewley grids

    hiimodel1 = interpol(hiigrids[*,logu1].nii_6584_oii,log12ohgrid,log12ohplot,/quad)
    hiimodel2 = interpol(hiigrids[*,logu2].nii_6584_oii,log12ohgrid,log12ohplot,/quad)
    sbmodel1 = interpol(sbgrids[*,logu1].nii_6584_oii,log12ohgrid,log12ohplot,/quad)
    sbmodel2 = interpol(sbgrids[*,logu2].nii_6584_oii,log12ohgrid,log12ohplot,/quad)    

    djs_oplot, log12ohplot, hiimodel1, line=0, thick=postthick, color=djs_icolor(talkcolor)
    djs_oplot, log12ohplot, hiimodel2, line=0, thick=postthick, color=djs_icolor(talkcolor)
;   djs_oplot, log12ohplot, sbmodel1, line=0, thick=postthick, color=djs_icolor(talkcolor)
;   djs_oplot, log12ohplot, sbmodel2, line=0, thick=postthick, color=djs_icolor(talkcolor)
    
    xyouts, log12ohplot[logu1], hiimodel1[logu1], string(U[logu1],format='(F5.2)'), /data, charsize=1.1, $
      charthick=postthick, align=1.0, color=djs_icolor(talkcolor)
    xyouts, log12ohplot[logu1], hiimodel2[logu1], string(U[logu2],format='(F5.2)'), /data, charsize=1.1, $
      charthick=postthick, align=1.0, color=djs_icolor(talkcolor)

;   legend, ['HII Region','Starburst'], line=[2,0], /right, /bottom, box=0, $
;     charsize=1.5, charthick=postthick, thick=postthick
       
    legend, '(a)', /left, /top, box=0, charsize=charsize_5, charthick=postthick, $
      textcolor=djs_icolor(talkcolor)
    
; --------------------------------------------------    
; Panel 2 - ([O III]/Hb)/([N II]/Ha)
; --------------------------------------------------    

    indx = where((hii.zstrong_oiiinii gt -900.0) and (hii.zt_log12oh_te gt -900.0),nindx)
    
    x = hii[indx].zt_log12oh_te
    xerr = hii[indx].zt_log12oh_te_err

    y = hii[indx].zstrong_oiiinii
    yerr = hii[indx].zstrong_oiiinii_err

    yrange = o3n2range2
    ytitle = 'log {([O III]/H\beta)/([N II]/H\alpha)}'
    
    zindicators_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, talk=dotalk, $
      xtitle='', ytitle='', ysty=11, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=3, /right, /top, position=pos[*,1], /noerase, ytickname=replicate(' ',10), $
      xtickname=replicate(' ',10), charsize=charsize_5, hiipsize=psize, /errorleft
    axis, /yaxis, ythick=postthick, ysty=1, yrange=!y.crange, ytitle=textoidl(ytitle), $
      charthick=postthick, charsize=charsize_5, color=djs_icolor(talkcolor)

; overplot the Pettini & Pagel (2004) and the Moustakas et
; al. relationships; only applicable from -1.0 to 1.9

;   o3n2_pettini = findgen((1.9-(-1.0))/0.01+1)*0.01+(-1.0)
;   oh_pettini = -0.32*o3n2_pettini + 8.73
;
;   o3n2_moustakas = findgen((1.9-(-1.0))/0.01+1)*0.01+(-1.0)
;   oh_moustakas = -0.2448*o3n2_moustakas + 8.5729
;
;   djs_oplot, oh_pettini, o3n2_pettini, line=2, thick=postthick+2, color='magenta'
;   djs_oplot, oh_moustakas, o3n2_moustakas, line=0, thick=postthick, color='green'

; overlay the kewley grids

    hiimodel1 = interpol(hiigrids[*,logu1].oiii_5007_h_beta_nii_6584_h_alpha,log12ohgrid,log12ohplot,/quad)
    hiimodel2 = interpol(hiigrids[*,logu2].oiii_5007_h_beta_nii_6584_h_alpha,log12ohgrid,log12ohplot,/quad)
    sbmodel1 = interpol(sbgrids[*,logu1].oiii_5007_h_beta_nii_6584_h_alpha,log12ohgrid,log12ohplot,/quad)
    sbmodel2 = interpol(sbgrids[*,logu2].oiii_5007_h_beta_nii_6584_h_alpha,log12ohgrid,log12ohplot,/quad)
    
    djs_oplot, log12ohplot, hiimodel1, line=0, thick=postthick, color=djs_icolor(talkcolor)
    djs_oplot, log12ohplot, hiimodel2, line=0, thick=postthick, color=djs_icolor(talkcolor)
;   djs_oplot, log12ohplot, sbmodel1, line=0, thick=postthick
;   djs_oplot, log12ohplot, sbmodel2, line=0, thick=postthick
    
    xyouts, log12ohplot[nlog12oh-1], hiimodel1[nlog12oh-1], string(U[logu1],format='(F5.2)'), /data, charsize=1.1, $
      charthick=postthick, align=0.0, color=djs_icolor(talkcolor)
    xyouts, log12ohplot[nlog12oh-1], hiimodel2[nlog12oh-1], string(U[logu2],format='(F5.2)'), /data, charsize=1.1, $
      charthick=postthick, align=0.0, color=djs_icolor(talkcolor)
;   xyouts, log12ohplot[logu1], hiimodel1[logu1], string(U[logu1],format='(F5.2)'), /data, charsize=1.1, charthick=postthick, align=1.0
;   xyouts, log12ohplot[logu1], hiimodel2[logu1], string(U[logu2],format='(F5.2)'), /data, charsize=1.1, charthick=postthick, align=1.0

;   legend, ['HII Region','Starburst'], line=[2,0], /left, /bottom, box=0, $
;     charsize=1.5, charthick=postthick, thick=postthick

    legend, '(b)', /left, /top, box=0, charsize=charsize_5, charthick=postthick, $
      textcolor=djs_icolor(talkcolor)

; --------------------------------------------------    
; Panel 3 - [N II]/Ha
; --------------------------------------------------    
 
    indx = where((hii.zstrong_niiha gt -900.0) and (hii.zt_log12oh_te gt -900.0),nindx)
    
    x = hii[indx].zt_log12oh_te
    xerr = hii[indx].zt_log12oh_te_err

    y = hii[indx].zstrong_niiha
    yerr = hii[indx].zstrong_niiha_err

    yrange = niiharange ; [-3,0]
    ytitle = 'log ([N II]/H\alpha)'
    
    zindicators_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, talk=dotalk, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=3, ystyle=3, /right, /top, position=pos[*,2], /noerase, $
      charsize=charsize_5, hiipsize=psize
 
; overplot the Pettini & Pagel (2004) and the Denicolo et al. (2002)
; relationships; only applicable from -2.5 to -0.3 in niiha; also
; overplot my best fit

;   n2_pettini = findgen(((-0.3)-(-2.5))/0.01+1)*0.01+(-2.5)
;   oh_pettini = 0.57*n2_pettini + 8.90
;
;   n2_denicolo = findgen((-0.3-(-2.5))/0.01+1)*0.01+(-2.5)
;   oh_denicolo = 0.73*n2_denicolo + 9.12
;
;   n2_moustakas = findgen((-0.5-(-2.5))/0.01+1)*0.01+(-2.5)
;   oh_moustakas = 0.794*n2_moustakas + 9.193
;
;   djs_oplot, oh_pettini, n2_pettini, line=2, thick=postthick+2, color='magenta'
;   djs_oplot, oh_denicolo, n2_denicolo, line=3, thick=postthick, color='red'
;   djs_oplot, oh_moustakas, n2_moustakas, line=0, thick=postthick, color='green'

; overlay some representative kewley grids

    hiimodel1 = interpol(hiigrids[*,logu1].nii_6584_h_alpha,log12ohgrid,log12ohplot,/quad)
    hiimodel2 = interpol(hiigrids[*,logu2].nii_6584_h_alpha,log12ohgrid,log12ohplot,/quad)
    sbmodel1 = interpol(sbgrids[*,logu1].nii_6584_h_alpha,log12ohgrid,log12ohplot,/quad)
    sbmodel2 = interpol(sbgrids[*,logu2].nii_6584_h_alpha,log12ohgrid,log12ohplot,/quad)
    
    djs_oplot, log12ohplot, hiimodel1, line=0, thick=postthick, color=djs_icolor(talkcolor)
    djs_oplot, log12ohplot, hiimodel2, line=0, thick=postthick, color=djs_icolor(talkcolor)
;   djs_oplot, log12ohplot, sbmodel1, line=0, thick=postthick, color=djs_icolor(talkcolor)
;   djs_oplot, log12ohplot, sbmodel2, line=0, thick=postthick, color=djs_icolor(talkcolor)
    
    xyouts, log12ohplot[nlog12oh-1], hiimodel1[nlog12oh-1], string(U[logu1],format='(F5.2)'), /data, charsize=1.1, $
      charthick=postthick, align=0.0, color=djs_icolor(talkcolor)
    xyouts, log12ohplot[nlog12oh-1], hiimodel2[nlog12oh-1], string(U[logu2],format='(F5.2)'), /data, charsize=1.1, $
      charthick=postthick, align=0.0, color=djs_icolor(talkcolor)
;   xyouts, log12ohplot[logu1], hiimodel1[logu1], string(U[logu1],format='(F5.2)'), /data, charsize=1.1, charthick=postthick, align=1.0
;   xyouts, log12ohplot[logu1], hiimodel2[logu1], string(U[logu2],format='(F5.2)'), /data, charsize=1.1, charthick=postthick, align=1.0
   
;   legend, ['HII Region','Starburst'], line=[2,0], /right, /bottom, box=0, $
;     charsize=1.5, charthick=postthick, thick=postthick

    legend, '(c)', /left, /top, box=0, charsize=charsize_5, charthick=postthick, $
      textcolor=djs_icolor(talkcolor)

; --------------------------------------------------    
; Panel 4 - [N II]/[S II]
; --------------------------------------------------    

    indx = where((hii.nii_6584_sii gt -900.0) and (hii.zt_log12oh_te gt -900.0),nindx)
    
    x = hii[indx].zt_log12oh_te
    xerr = hii[indx].zt_log12oh_te_err

    y = hii[indx].nii_6584_sii
    yerr = hii[indx].nii_6584_sii_err

    yrange = niisiirange
    ytitle = 'log ([N II]/[S II])'
    
    zindicators_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, talk=dotalk, $
      xtitle=xtitle, ytitle='', ystyle=11, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=3, /right, /top, position=pos[*,3], /noerase, $
      ytickname=replicate(' ',10), charsize=charsize_5, hiipsize=psize, /errorleft
    axis, /yaxis, ythick=postthick, ysty=1, yrange=!y.crange, ytitle=textoidl(ytitle), $
      charthick=postthick, charsize=charsize_5, color=djs_icolor(talkcolor)

; overlay some representative kewley grids

    hiimodel1 = interpol(hiigrids[*,logu1].nii_6584_sii,log12ohgrid,log12ohplot,/quad)
    hiimodel2 = interpol(hiigrids[*,logu2].nii_6584_sii,log12ohgrid,log12ohplot,/quad)
    sbmodel1 = interpol(sbgrids[*,logu1].nii_6584_sii,log12ohgrid,log12ohplot,/quad)
    sbmodel2 = interpol(sbgrids[*,logu2].nii_6584_sii,log12ohgrid,log12ohplot,/quad)
    
    djs_oplot, log12ohplot, hiimodel1, line=0, thick=postthick, color=djs_icolor(talkcolor)
    djs_oplot, log12ohplot, hiimodel2, line=0, thick=postthick, color=djs_icolor(talkcolor)
;   djs_oplot, log12ohplot, sbmodel1, line=0, thick=postthick, color=djs_icolor(talkcolor)
;   djs_oplot, log12ohplot, sbmodel2, line=0, thick=postthick, color=djs_icolor(talkcolor)
    
    xyouts, log12ohplot[nlog12oh-1], hiimodel1[nlog12oh-1], string(U[logu1],format='(F5.2)'), /data, charsize=1.1, $
      charthick=postthick, align=0.0, color=djs_icolor(talkcolor)
    xyouts, log12ohplot[nlog12oh-1], hiimodel2[nlog12oh-1], string(U[logu2],format='(F5.2)'), /data, charsize=1.1, $
      charthick=postthick, align=0.0, color=djs_icolor(talkcolor)
;   xyouts, log12ohplot[logu1], hiimodel1[logu1], string(U[logu1],format='(F5.2)'), /data, charsize=1.1, charthick=postthick, align=1.0
;   xyouts, log12ohplot[logu1], hiimodel2[logu1], string(U[logu2],format='(F5.2)'), /data, charsize=1.1, charthick=postthick, align=1.0

;   legend, ['HII Region','Starburst'], line=[2,0], /right, /top, box=0, $
;     charsize=1.5, charthick=postthick, thick=postthick

    legend, '(d)', /left, /top, box=0, charsize=charsize_5, charthick=postthick, $
      textcolor=djs_icolor(talkcolor)

    im_openclose, postscript=(keyword_set(postscript) or keyword_set(pdf)), /close    
    if keyword_set(pdf) then begin
       splog, 'Writing '+pspath+psname+'.pdf'
       spawn, 'ps2pdf13 '+pspath+psname+'.ps '+pspath+psname+'.pdf', /sh
       rmfile, pspath+psname+'.ps' ; delete the junk ps file
    endif

; ------------------------------------------------------------
; Te [direct] vs Te [literature] 
; ------------------------------------------------------------

    psname = 'hii_te_me_vs_te_lit'
    im_openclose, pspath+psname, postscript=postscript, encapsulated=encapsulated, xsize=8.5, ysize=8.0

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.6, height=6.6, $
      xmargin=[1.5,0.4], ymargin=[0.3,1.1], xpage=8.5, ypage=8.0, $
      position=pos, /normal

    indx = where((hii.lit_t4363 gt -900.0) and (hii.zt_t4363 gt -900.0),nindx)
;   indx = where((hii.lit_t4363 gt -900.0) and (hii.zt_t4363 gt -900.0) and $
;     (hii.lit_t4363-hii.zt_t4363 ne 0.0),nindx)
    
    x = hii[indx].zt_t4363/1E4
    xerr = hii[indx].zt_t4363_err/1E4

    y = hii[indx].lit_t4363/1E4
    yerr = hii[indx].lit_t4363_err/1E4

    stats = im_stats((y-x)*1E4,/verbose,/baremin)

    xrange = Terange/1E4
    yrange = xrange

    xtitle = 'Electron Temperature [10^{4} K]'
    ytitle = 'Electron Temperature [Literature, 10^{4} K]'

    zindicators_lineplot, x, y, xerr, yerr, plottype=2, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      hiipsize=1.0, position=pos[*,0], charsize=singlecharsize
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick, color=talkcolor
    legend, '(a)', /left, /top, box=0, charsize=singlecharsize, charthick=postthick
    
    im_openclose, postscript=(keyword_set(postscript) or keyword_set(pdf)), /close    
    if keyword_set(pdf) then begin
       splog, 'Writing '+pspath+psname+'.pdf'
       spawn, 'ps2pdf13 '+pspath+psname+'.ps '+pspath+psname+'.pdf', /sh
       rmfile, pspath+psname+'.ps' ; delete the junk ps file
    endif

; ------------------------------------------------------------
; 12+log(O/H) [Me] vs 12+log(O/H) [Literature]
; ------------------------------------------------------------

    psname = 'hii_12oh_me_vs_12oh_lit'
    im_openclose, pspath+psname, postscript=postscript, encapsulated=encapsulated, xsize=8.5, ysize=8.0

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.6, height=6.6, $
      xmargin=[1.5,0.4], ymargin=[0.3,1.1], xpage=8.5, ypage=8.0, $
      position=pos, /normal

    indx = where((hii.lit_log12oh_te gt -900.0) and (hii.zt_log12oh_te gt -900.0),nindx)
    
    x = hii[indx].zt_log12oh_te
    xerr = hii[indx].zt_log12oh_te_err

    y = hii[indx].lit_log12oh_te
    yerr = hii[indx].lit_log12oh_te_err

    stats = im_stats(y-x,/verbose,/baremin)
    
    xrange = ohrange2
    yrange = xrange

    xtitle = '12 + log (O/H)'
    ytitle = '12 + log (O/H) [Literature]'

    zindicators_lineplot, x, y, xerr, yerr, plottype=2, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      hiipsize=1.0, position=pos[*,0], charsize=singlecharsize
    djs_oplot, !x.crange, !y.crange, line=0, thick=postthick, color=talkcolor
    legend, '(b)', /left, /top, box=0, charsize=singlecharsize, charthick=postthick

    im_openclose, postscript=(keyword_set(postscript) or keyword_set(pdf)), /close    
    if keyword_set(pdf) then begin
       splog, 'Writing '+pspath+psname+'.pdf'
       spawn, 'ps2pdf13 '+pspath+psname+'.ps '+pspath+psname+'.pdf', /sh
       rmfile, pspath+psname+'.ps' ; delete the junk ps file
    endif

; --------------------------------------------------    
; SELECT PLOTS FOR THE PAPER HERE
; --------------------------------------------------    

    if keyword_set(paper) then begin

       splog, 'Writing paper plots to '+paperpath+'.'
       paperplots = [$
         'hii_ratios_vs_12oh_te_3panel',$
         'hii_12oh_te_vs_oiiiratios_4panel',$
         'hii_12oh_te_vs_niiratios_4panel',$
         'hii_te_me_vs_te_lit',$
         'hii_12oh_me_vs_12oh_lit'$
         ]+'.eps'

       for k = 0L, n_elements(paperplots)-1L do $
         spawn, ['/bin/cp -f '+html_path+htmlbase+'/'+paperplots[k]+' '+paperpath], /sh

    endif

; --------------------------------------------------    
; GENERATE THE WEB PAGE
; --------------------------------------------------    
    
    if keyword_set(postscript) and keyword_set(encapsulated) then $
      im_ps2html, htmlbase, html_path=html_path, cleanpng=0, npscols=3, _extra=extra
    
stop    
    
; ------------------------------------------------------------
; ([O III]/Hb)/([N II]/Ha) vs 12+log (O/H) [Te] - Pettini & Pagel (2004)
; ------------------------------------------------------------

    psname = 'hii_oiiinii_vs_12oh_te'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.0, /encapsulated

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.6, height=6.6, $
      xmargin=[1.5,0.4], ymargin=[0.3,1.1], xpage=8.5, ypage=8.0, $
      position=pos, /normal

    indx = where((hii.oiii_5007_h_beta_nii_6584_h_alpha gt -900.0) and (hii.zt_log12oh_te gt -900.0),nindx)
    
    x = hii[indx].oiii_5007_h_beta_nii_6584_h_alpha
    xerr = hii[indx].oiii_5007_h_beta_nii_6584_h_alpha_err

    y = hii[indx].zt_log12oh_te
    yerr = hii[indx].zt_log12oh_te_err

    xrange = o3n2range
    yrange = ohrange
    
;   xtitle = 'log {([O III]/H\beta)/([N II]/H\alpha)}'
    xtitle = 'log {([O III] \lambda5007/H\beta)/([N II] \lambda6584/H\alpha)}'
    ytitle = '12 + log (O/H) [T_{e}]'
    
    zindicators_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      hiipsize=1.0, position=pos[*,0], charsize=charsize_8

; fit a linear function/bisector

    range = where((x ge min(x)) and (x le 1.9))
    fitx = x[range]
    fity = y[range]
    fitxerr = xerr[range]
    fityerr = yerr[range]

    coeff = poly_fit(fitx,fity,1,sigma=coeff_err,yfit=smallfit);,measure_errors=fityerr)

;   sixlin, fitx, fity, a, siga, b, sigb
;   coeff = [a[2],b[2]]
;   coeff_err = [siga[2],sigb[2]]

; note I changed min(fitx) to -1.0!

    xaxis = findgen(((1.9)-(-1.0))/0.001+1)*0.001+(-1.0)
;   xaxis = findgen((max(fitx)-min(fitx))/0.001+1)*0.001+min(fitx)
    yfit = poly(xaxis,coeff)
    djs_oplot, xaxis, yfit, line=0, thick=postthick, color='green'

    splog, 'oiiinii vs 12+log(O/H) [Te] coefficients and range: '
    niceprint, coeff, coeff_err, minmax(fitx)

; overplot the Pettini & Pagel (2004) relationship; only applicable
; from -1.0 to 1.9 in oiiinii

    o3n2_pettini = findgen((1.9-(-1.0))/0.01+1)*0.01+(-1.0)
    oh_pettini = -0.32*o3n2_pettini + 8.73 
    djs_oplot, o3n2_pettini, oh_pettini, line=2, thick=postthick+2, color=fsc_color('orchid',189)

    legend, ['This Paper','Pettini & Pagel 2004'], line=[0,2], thick=postthick, /right, /top, box=0, $
      charsize=1.5, charthick=postthick, color=djs_icolor(['green','magenta']), $
      textcolor=djs_icolor(talkcolor)

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; 12+log (O/H) [Te] vs [O III]/Hb
; ------------------------------------------------------------

    psname = 'hii_oiiihb_vs_12oh_te'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.0, /encapsulated

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.6, height=6.6, $
      xmargin=[1.5,0.4], ymargin=[0.3,1.1], xpage=8.5, ypage=8.0, $
      position=pos, /normal

    indx = where((hii.oiii_5007_h_beta gt -900.0) and (hii.zt_log12oh_te gt -900.0),nindx)
    
    x = hii[indx].oiii_5007_h_beta
    xerr = hii[indx].oiii_5007_h_beta_err

    y = hii[indx].zt_log12oh_te
    yerr = hii[indx].zt_log12oh_te_err

    xrange = oiiihbrange
    yrange = ohrange

    xtitle = 'log ([O III] \lambda5007 / H\beta)'
    ytitle = '12 + log (O/H) [T_{e}]'

    zindicators_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      hiipsize=1.0, position=pos[*,0], charsize=charsize_8

; fit a linear function/bisector

    range = where((x ge -1.5) and (x le 0.3))
    fitx = x[range]
    fity = y[range]
    fitxerr = xerr[range]
    fityerr = yerr[range]

    coeff = poly_fit(fitx,fity,1,sigma=coeff_err);,measure_errors=fityerr)

;   sixlin, fitx, fity, a, siga, b, sigb
;   coeff = [a[2],b[2]]
;   coeff_err = [siga[2],sigb[2]]

    xaxis = findgen((max(fitx)-min(fitx))/0.001+1)*0.001+min(fitx)
    yfit = poly(xaxis,coeff)
;   djs_oplot, xaxis, yfit, line=0, thick=postthick, color='dark green'

;   splog, '[O III]/Hb vs 12+log(O/H) [Te] coefficients: '
;   niceprint, coeff, coeff_err, minmax(fitx)

; overplot the Melbourne & Salzer relation; only applicable from -1.0
; to 0.7 in [O III]/Hb

    oiiihb_melbourne = findgen((0.7-(-1.0))/0.01+1)*0.01+(-1.0)
    oh_melbourne = -0.663*oiiihb_melbourne + 8.65
    djs_oplot, oiiihb_melbourne, oh_melbourne, line=5, thick=postthick, color='orange'

    legend, ['Melbourne & Salzer 2002'], line=5, thick=postthick, /right, /top, box=0, $
      charsize=1.5, charthick=postthick, color=djs_icolor('orange'), $
      textcolor=djs_icolor(talkcolor)
;   legend, ['This Paper','Melbourne & Salzer 2002'], line=[0,5], thick=postthick, /right, /top, box=0, $
;     charsize=1.5, charthick=postthick, color=djs_icolor(['dark green','orange']), $
;     textcolor=djs_icolor(talkcolor)

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; [N II]/Ha vs 12+log (O/H) [Te] 
; ------------------------------------------------------------

; van Zee et al. (1998), Denicolo et al. (2002), Storchi-Bergmann,
; Calzetti, & Kinney (1994), Raimann et al. (2000), Pettini & Pagel
; (2004) 
    
    psname = 'hii_niiha_vs_12oh_te'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.0, /encapsulated

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.6, height=6.6, $
      xmargin=[1.5,0.4], ymargin=[0.3,1.1], xpage=8.5, ypage=8.0, $
      position=pos, /normal

    indx = where((hii.zstrong_niiha gt -900.0) and (hii.zt_log12oh_te gt -900.0),nindx)
    
    x = hii[indx].zstrong_niiha
    xerr = hii[indx].zstrong_niiha_err

    y = hii[indx].zt_log12oh_te
    yerr = hii[indx].zt_log12oh_te_err

    xrange = niiharange
    yrange = ohrange

    xtitle = 'log ([N II] \lambda6584 / H\alpha)'
    ytitle = '12 + log (O/H) [T_{e}]'
    
    zindicators_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      hiipsize=1.0, position=pos[*,0], charsize=charsize_8

; overplot various calibrations

    niiha_denicolo = findgen((0.0-(-2.6))/0.01+1)*0.01+(-2.6)
    ohaxis_denicolo = 0.73*niiha_denicolo + 9.12

    niiha_melbourne = findgen((-0.45-min(x))/0.01+1)*0.01+min(x)
    ohaxis_melbourne = 9.26 + 1.23*niiha_melbourne + 0.204*niiha_melbourne^2
    
;   ohaxis_vanzee = findgen((9.1-7.0)/0.01+1)*0.01+7.0
;   niiha_vanzee = (ohaxis_vanzee-9.36) / 1.02

    ohaxis_pettini = findgen((8.729-7.475)/0.01+1)*0.01+7.475 ; only valid in this range!
    niiha_pettini = (ohaxis_pettini-8.90) / 0.57

    djs_oplot, niiha_denicolo, ohaxis_denicolo, line=3, thick=postthick, color='red'
    djs_oplot, niiha_pettini, ohaxis_pettini, line=2, thick=postthick+2, color='magenta'
;   djs_oplot, niiha_melbourne, ohaxis_melbourne, line=5, thick=postthick, color='orange'

    legend, ['This Paper','Denicolo et al. 2002','Pettini & Pagel 2004'], $
      line=[0,3,2], thick=postthick, $
      /left, /top, box=0, charsize=1.5, charthick=postthick, $
      color=djs_icolor(['green','red','magenta']), $
      textcolor=djs_icolor(talkcolor)

;   legend, ['This Paper','Denicolo et al. 2002','Melbourne & Salzer 2002','Pettini & Pagel 2004'], $
;     line=[0,3,5,2], thick=postthick, $
;     /left, /top, box=0, charsize=1.5, charthick=postthick, $
;     color=djs_icolor(['dark green','red','orange','purple']), $
;     textcolor=djs_icolor(talkcolor)

; fit a linear function/bisector 

    niihamin = min(x)
    niihamax = -0.5
    
    range = where((x ge niihamin) and (x le niihamax))
    fitx = x[range]
    fity = y[range]
    fitxerr = xerr[range]
    fityerr = yerr[range]

;   coeff = poly_fit(fitx,fity,1,sigma=coeff_err);,measure_errors=fityerr)

    sixlin, fitx, fity, a, siga, b, sigb
    coeff = [a[2],b[2]]
    coeff_err = [siga[2],sigb[2]]

    xaxis = findgen((max(fitx)-min(fitx))/0.001+1)*0.001+min(fitx)
    yfit = poly(xaxis,coeff)
    djs_oplot, xaxis, yfit, line=0, thick=postthick, color='green'

    splog, '[N II]/Ha vs 12+log(O/H) [Te] coefficients: '
    niceprint, coeff, coeff_err, minmax(fitx)

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; O32 vs 12+log (O/H) [Te]
; ------------------------------------------------------------

    psname = 'hii_o32_vs_12oh_te'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.0, /encapsulated

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.6, height=6.6, $
      xmargin=[1.5,0.4], ymargin=[0.3,1.1], xpage=8.5, ypage=8.0, $
      position=pos, /normal

    indx = where((hii.zstrong_o32 gt -900.0) and (hii.zt_log12oh_te gt -900.0),nindx)
    
    x = hii[indx].zstrong_o32
    xerr = hii[indx].zstrong_o32_err

    y = hii[indx].zt_log12oh_te
    yerr = hii[indx].zt_log12oh_te_err

    xrange = o32range
    yrange = ohrange

    xtitle = 'log (O_{32})'
;   xtitle = 'log ([O III] \lambda\lambda4959,5007 / [O II] \lambda3727)'
    ytitle = '12 + log (O/H) [T_{e}]'
    
    zindicators_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      hiipsize=1.0, position=pos[*,0], charsize=charsize_8

; fit a linear function/bisector

    range = where((x ge min(x)) and (x le max(x)))
    fitx = x[range]
    fity = y[range]
    fitxerr = xerr[range]
    fityerr = yerr[range]

    coeff = poly_fit(fitx,fity,1,sigma=coeff_err);,measure_errors=fityerr)

;   sixlin, fitx, fity, a, siga, b, sigb
;   coeff = [a[2],b[2]]
;   coeff_err = [siga[2],sigb[2]]

    xaxis = findgen((max(fitx)-min(fitx))/0.001+1)*0.001+min(fitx)
    yfit = poly(xaxis,coeff)
;   djs_oplot, xaxis, yfit, line=0, thick=postthick, color='dark green'
;
;   splog, 'o32 vs 12+log(O/H) [Te] coefficients: '
;   niceprint, coeff, coeff_err, minmax(fitx)
;      
;   legend, 'This Paper', line=0, thick=postthick, /right, /top, box=0, $
;     charsize=1.5, charthick=postthick, color=djs_icolor('dark green')

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; P vs 12+log (O/H) [Te]
; ------------------------------------------------------------

    psname = 'hii_P_vs_12oh_te'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.0, /encapsulated

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.6, height=6.6, $
      xmargin=[1.5,0.4], ymargin=[0.3,1.1], xpage=8.5, ypage=8.0, $
      position=pos, /normal

    indx = where((hii.zstrong_P gt -900.0) and (hii.zt_log12oh_te gt -900.0),nindx)
    
    x = hii[indx].zstrong_P
    xerr = hii[indx].zstrong_P_err

    y = hii[indx].zt_log12oh_te
    yerr = hii[indx].zt_log12oh_te_err

    xrange = prange
    yrange = ohrange

    xtitle = 'P'
;   xtitle = 'log {[O II] \lambda3727 / ([O II] \lambda3727 + [O III] \lambda\lambda4959,5007)}'
    ytitle = '12 + log (O/H) [T_{e}]'
    
    zindicators_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      hiipsize=1.0, position=pos[*,0], charsize=charsize_8

; fit a linear function/bisector

    range = where((x ge min(x)) and (x le max(x)))
    fitx = x[range]
    fity = y[range]
    fitxerr = xerr[range]
    fityerr = yerr[range]

    coeff = poly_fit(fitx,fity,1,sigma=coeff_err);,measure_errors=fityerr)

;   sixlin, fitx, fity, a, siga, b, sigb
;   coeff = [a[2],b[2]]
;   coeff_err = [siga[2],sigb[2]]

;   medxbin, fitx, fity, 0.1, medx=medx, medy=medy, sigy=sigy, thresh=0

    xaxis = findgen((max(fitx)-min(fitx))/0.001+1)*0.001+min(fitx)
    yfit = poly(xaxis,coeff)
;   djs_oplot, xaxis, yfit, line=0, thick=postthick, color='dark green'
;
;   splog, 'P vs 12+log(O/H) [Te] coefficients: '
;   niceprint, coeff, coeff_err, minmax(fitx)
;      
;   legend, 'This Paper', line=0, thick=postthick, /right, /top, box=0, $
;     charsize=1.5, charthick=postthick, color=djs_icolor('dark green')

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; [N II]/[S II] vs 12+log (O/H) [Te]
; ------------------------------------------------------------

    psname = 'hii_niisii_vs_12oh_te'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.0, /encapsulated

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.6, height=6.6, $
      xmargin=[1.5,0.4], ymargin=[0.3,1.1], xpage=8.5, ypage=8.0, $
      position=pos, /normal

    indx = where((hii.nii_6584_sii gt -900.0) and (hii.zt_log12oh_te gt -900.0),nindx)
    
    x = hii[indx].nii_6584_sii
    xerr = hii[indx].nii_6584_sii_err

    y = hii[indx].zt_log12oh_te
    yerr = hii[indx].zt_log12oh_te_err

    xrange = niisiirange
    yrange = ohrange

    xtitle = 'log ([N II] \lambda6584 / [S II] \lambda\lambda6716,31)'
    ytitle = '12 + log (O/H) [T_{e}]'
    
    zindicators_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      hiipsize=1.0, position=pos[*,0], charsize=charsize_8

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; [N II]/[O II] vs 12+log (O/H) [Te]
; ------------------------------------------------------------

    psname = 'hii_niioii_vs_12oh_te'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.0, /encapsulated

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.6, height=6.6, $
      xmargin=[1.5,0.4], ymargin=[0.3,1.1], xpage=8.5, ypage=8.0, $
      position=pos, /normal

    indx = where((hii.zstrong_niioii gt -900.0) and (hii.zt_log12oh_te gt -900.0),nindx)
    
    x = hii[indx].zstrong_niioii
    xerr = hii[indx].zstrong_niioii_err

    y = hii[indx].zt_log12oh_te
    yerr = hii[indx].zt_log12oh_te_err

    xrange = niioiirange
    yrange = ohrange

    xtitle = 'log ([N II] \lambda6584 / [O II])'
    ytitle = '12 + log (O/H) [T_{e}]'
    
    zindicators_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      hiipsize=1.0, position=pos[*,0], charsize=charsize_8

; fit a linear function/bisector

    range = where((x ge -1.0) and (x le max(x)))
    fitx = x[range]
    fity = y[range]
    fitxerr = xerr[range]
    fityerr = yerr[range]

    coeff = poly_fit(fitx,fity,1,sigma=coeff_err);,measure_errors=fityerr)

;   sixlin, fitx, fity, a, siga, b, sigb
;   coeff = [a[2],b[2]]
;   coeff_err = [siga[2],sigb[2]]

    xaxis = findgen((max(fitx)-min(fitx))/0.001+1)*0.001+min(fitx)
    yfit = poly(xaxis,coeff)
    djs_oplot, xaxis, yfit, line=0, thick=postthick, color='green'

    splog, '[N II]/[O II] vs 12+log(O/H) [Te] coefficients: '
    niceprint, coeff, coeff_err, minmax(fitx)

; overplot the Kewley & Dopita (2002) relation for high metallicity
    
    niioii_kd02 = findgen((1.0-(-0.97))/0.01+1)*0.01-0.97
    ohaxis_kd02 = alog10(1.54020 + 1.26602*niioii_kd02 + 0.167977*niioii_kd02^2) + 8.93 - 0.24 ; <-- NOTE!!!
    djs_oplot, niioii_kd02, ohaxis_kd02, line=6, thick=postthick, color='cyan'
    
    legend, ['This Paper','Kewley & Dopita 2002'], line=[0,6], thick=postthick, /left, /top, box=0, $
      charsize=1.5, charthick=postthick, color=djs_icolor(['green','cyan']), $
      textcolor=djs_icolor(talkcolor)

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; R23 vs 12+log (O/H) [Te]
; ------------------------------------------------------------

    psname = 'hii_r23_vs_12oh_te'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.0, /encapsulated

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.6, height=6.6, $
      xmargin=[1.5,0.4], ymargin=[0.3,1.1], xpage=8.5, ypage=8.0, $
      position=pos, /normal

    indx = where((hii.zstrong_r23 gt -900.0) and (hii.zt_log12oh_te gt -900.0),nindx)
     
    x = alog10(hii[indx].zstrong_r23)
    xerr = hii[indx].zstrong_r23_err/hii[indx].zstrong_r23/alog(10.0)

    y = hii[indx].zt_log12oh_te
    yerr = hii[indx].zt_log12oh_te_err

    xrange = r23range
    yrange = ohrange

    xtitle = 'log (R_{23})'
    ytitle = '12 + log (O/H) [T_{e}]'
    
    zindicators_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      hiipsize=1.0, position=pos[*,0], charsize=charsize_8

; overplot the Kewley and McGaugh theoretical calibrations
 
    logo32 = 0.0
    logq = alog10(4E7)
    logr23 = findgen((1.0-(-0.5))/0.01+1)*0.01-0.5
 
    M91 = 12.0 - 2.939 - 0.2*logr23 - 0.237*logr23^2 - 0.305*logr23^3 - $
      0.0283*logr23^4 - logo32*(0.0047 - 0.0221*logr23 - 0.102*logr23^2 - $
      0.0817*logr23^3 - 0.00717*logr23^4)
 
    KK04 = 9.72D - 0.777*logr23 - 0.951*logr23^2 - 0.072*logr23^3 - 0.811*logr23^4 - $
      logq*(0.0737 - 0.0713*logr23 - 0.141*logr23^2 + 0.0373*logr23^3 - 0.058*logr23^4)
 
    djs_oplot, logr23, M91, line=0, thick=postthick, color='red'
    djs_oplot, logr23, KK04, line=2, thick=postthick, color='green'

    legend, ['McGaugh (1991)','Kobulnicky & Kewley (2004)'], $
      line=[0,2], thick=postthick, /left, /bottom, box=0, charsize=1.5, $
      charthick=postthick, color=djs_icolor(['red','green'])
    
; overplot the Pilyugin (2000,2001) parameterizations

;; predict R3 based on a robust linefit to the hii regions    
;
;   r23max = 1.0 & r23min = 0.3
;   
;   good = where((hii.zstrong_r23 gt r23min) and (hii.zstrong_r23 lt r23max) and (hii.zstrong_R3 gt -900))
;   r23 = 10^hii[good].zstrong_r23
;   R3 = 10^hii[good].zstrong_R3
;   coeff = robust_poly_fit(r23,R3,2,yfit)
;   srt = sort(r23)
;   plot, r23, R3, ps=4, xsty=3, ysty=3
;   djs_oplot, r23[srt], yfit[srt], line=0, thick=2
;
; lower branch    
;   
;   r23 = findgen((r23max-r23min)/0.01+1)*0.01+r23min
;   R3 = coeff[0] + coeff[1]*r23 + coeff[2]*r23^2
;   
;   oh = 6.35 + 3.19*r23 - 1.74*R3
;   keep = where(oh lt 7.95,nkeep)
;   
;   djs_oplot, oh[keep], r23[keep], line=1, thick=postthick, color='orange'
;
; upper branch    
;
;   P = 0.3
;   r23 = 10^(findgen((1.0-(-0.5))/0.01+1)*0.01-0.5)
;   oh = (r23 + 54.2 + 59.45*P + 7.31*P^2) / (6.07 + 6.71*P + 0.37*P^2 + 0.243*r23)
;   keep = where(oh gt 8.2,nkeep)
;   
;   djs_oplot, alog10(r23[keep]), oh[keep], line=5, thick=postthick, color='orange'
;
;   legend, ['McGaugh (1991)','Pilyugin & Thuan (2005, P=0.3)','Kobulnicky & Kewley (2004)'], $
;     line=[0,5,2], thick=postthick, /left, /bottom, box=0, charsize=1.5, $
;     charthick=postthick, color=djs_icolor(['red','orange','green']), $
;     textcolor=djs_icolor(talkcolor)
;   
;; fit a polynomial to the median distribution - upper branch
;;   
;;   range = where((x le 1.0) and (y gt 8.3))
;;   fitx = x[range]
;;   fity = y[range]
;;
;;   srt = sort(fitx)
;;   
;;   med = dkboxstats(fity[srt],boxstat='median',xwidth=5,lower=losig,upper=hisig)
;;   sigma = dkboxstats(fity[srt],boxstat='sigma',xwidth=5,sigrej=3.0)
;;
;;   plotsym, 8, 1.0
;;   oploterror, fitx[srt], med, sig, /nohat, color='red'
;;   djs_oplot, fitx[srt], med, color='red', line=0, thick=postthick
;;   djs_oplot, fitx[srt], med+hisig, color='red', line=2, thick=postthick
;;   djs_oplot, fitx[srt], med-losig, color='red', line=2, thick=postthick
;;   
;;   coeff = poly_fit(fitx[srt],med,2,measure_errors=sigma,sigma=coeff_err)
;;   xaxis = findgen((max(fitx[srt])-min(fitx[srt]))/0.01+1)*0.01+min(fitx[srt])
;;   medfit = poly(xaxis,coeff)
;;   djs_oplot, xaxis, medfit, line=0, thick=postthick+2, color='dark green'
;;   
;;   splog, 'r23 - upper - vs 12+log(O/H) [Te] coefficients: '
;;   niceprint, coeff, coeff_err

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; 12+log(O/H) Te vs 12+log(O/H) Strong
; ------------------------------------------------------------

    psname = 'hii_12oh_te_vs_12oh_strong_3panel'
    im_openclose, pspath+psname, postscript=postscript, encapsulated=encapsulated, xsize=8.5, ysize=10.5

    pagemaker, nx=1, ny=3, height=3.0*[1,1,1], width=6.8, xmargin=[1.3,0.4], $
      ymargin=[0.4,1.1], xspace=0, yspace=0, xpage=8.5, ypage=10.5, $
      position=pos, /normal

    if keyword_set(dotalk) then polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, $
      color=djs_icolor('black')

    ytitle = ohtitle
    xtitle = ohtitle+' T_{e}'

    xrange = ohrange7
    yrange = ohrange7

; --------------------------------------------------    
; Panel 1 - 12+log(O/H) Te vs 12+log(O/H) M91
; --------------------------------------------------    

    indx = where((hii.zt_log12oh_te gt -900.0) and (hii.zstrong_12oh_m91_upper gt -900.0) and $
      (hii.zstrong_12oh_m91_lower gt -900.0),nindx)

    x = hii[indx].zt_log12oh_te
    xerr = hii[indx].zt_log12oh_te_err

    y = x*0.0
    yerr= xerr*0.0
    
    for idiv = 0L, ndiv-1L do begin

       up = where((hii[indx].zt_log12oh_te gt div[idiv]),nup,comp=lo,ncomp=nlo)
       if (nup ne 0L) then begin
          y[up] = hii[indx[up]].zstrong_12oh_m91_upper
          yerr[up] = hii[indx[up]].zstrong_12oh_m91_upper_err
       endif

       if (nlo ne 0L) then begin
          y[lo] = hii[indx[lo]].zstrong_12oh_m91_lower
          yerr[lo] = hii[indx[lo]].zstrong_12oh_m91_lower_err
       endif

       resid[idiv] = stddev(x-y)
       
    endfor

    minresid = min(resid,minindx)
    mindiv = div[minindx]
    splog, '12+log(O/H) Te vs M91: '+string(mindiv,format='(F4.2)')

    up = where((hii[indx].zt_log12oh_te gt mindiv),nup)
    if (nup ne 0L) then begin
       xup = hii[indx[up]].zt_log12oh_te
       xerrup = hii[indx[up]].zt_log12oh_te_err

       yup = hii[indx[up]].zstrong_12oh_m91_upper
       yerrup = hii[indx[up]].zstrong_12oh_m91_upper_err

       x = xup & xerr = xerrup
       y = yup & yerr = yerrup
    endif

    lo = where((hii[indx].zt_log12oh_te lt mindiv),nlo)
    if (nlo ne 0L) then begin
       xlo = hii[indx[lo]].zt_log12oh_te
       xerrlo = hii[indx[lo]].zt_log12oh_te_err

       ylo = hii[indx[lo]].zstrong_12oh_m91_lower
       yerrlo = hii[indx[lo]].zstrong_12oh_m91_lower_err

       x = [xup,xlo] & xerr = [xerrup,xerrlo]
       y = [yup,ylo] & yerr = [yerrup,yerrlo]
    endif

    upstats = im_stats(xup-yup)
    upxstr = strtrim(string(upstats.median,format='(F12.2)'),2)+'\pm'+$
      strtrim(string(upstats.sig68mean,format='(F12.2)'),2)

    lostats = im_stats(xlo-ylo)
    loxstr = strtrim(string(lostats.median,format='(F12.2)'),2)+'\pm'+$
      strtrim(string(lostats.sig68mean,format='(F12.2)'),2)

    stats = im_stats(x-y,/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)

    zindicators_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, talk=dotalk, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, charsize=charsize_6, $
      xstyle=3, ystyle=3, /right, /top, position=pos[*,0], xtickname=replicate(' ',10), hiipsize=0.8
    djs_oplot, !x.crange, !y.crange, thick=postthick, color=djs_icolor(talkcolor), line=2
 
    legend, textoidl([upxstr+' (Upper)',loxstr+' (Lower)']), /right, /bottom, box=0, $
      charsize=charsize_5, charthick=postthick, textcolor=djs_icolor(talkcolor)
    legend, 'McGaugh 1991', /left, /top, box=0, charsize=charsize_5, charthick=postthick, $
      textcolor=djs_icolor(talkcolor)

; --------------------------------------------------    
; Panel 2 - 12+log(O/H) Te vs 12+log(O/H) ZKH94
; --------------------------------------------------    

    indx = where((hii.zt_log12oh_te gt 8.0) and (hii.zstrong_12oh_zkh94 gt -900.0),nindx)

    x = hii[indx].zt_log12oh_te
    xerr = hii[indx].zt_log12oh_te_err

    y = hii[indx].zstrong_12oh_zkh94
    yerr = hii[indx].zstrong_12oh_zkh94_err

    stats = im_stats(x-y,/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)

    zindicators_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, talk=dotalk, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, charsize=charsize_6, $
      xstyle=3, ystyle=3, /right, /top, position=pos[*,1], xtickname=replicate(' ',10), hiipsize=0.8, $
      /noerase
    djs_oplot, !x.crange, !y.crange, thick=postthick, color=djs_icolor(talkcolor), line=2
 
    legend, textoidl(xstr+' (Upper)'), /right, /bottom, box=0, $
      charsize=charsize_5, charthick=postthick, textcolor=djs_icolor(talkcolor)
    legend, 'Zaritsky et al. 1994', /left, /top, box=0, charsize=charsize_5, charthick=postthick, $
      textcolor=djs_icolor(talkcolor)

; --------------------------------------------------    
; Panel 3 - 12+log(O/H) Te vs 12+log(O/H) KK04
; --------------------------------------------------    

    indx = where((hii.zt_log12oh_te gt -900.0) and (hii.zstrong_12oh_kk04_r23_upper gt -900.0) and $
      (hii.zstrong_12oh_kk04_r23_lower gt -900.0),nindx)

    x = hii[indx].zt_log12oh_te
    xerr = hii[indx].zt_log12oh_te_err

    y = x*0.0
    yerr= xerr*0.0
    
    for idiv = 0L, ndiv-1L do begin

       up = where((hii[indx].zt_log12oh_te gt div[idiv]),nup,comp=lo,ncomp=nlo)
       if (nup ne 0L) then begin
          y[up] = hii[indx[up]].zstrong_12oh_kk04_r23_upper
          yerr[up] = hii[indx[up]].zstrong_12oh_kk04_r23_upper_err
       endif

       if (nlo ne 0L) then begin
          y[lo] = hii[indx[lo]].zstrong_12oh_kk04_r23_lower
          yerr[lo] = hii[indx[lo]].zstrong_12oh_kk04_r23_lower_err
       endif

       resid[idiv] = stddev(x-y)
       
    endfor

    minresid = min(resid,minindx)
    mindiv = div[minindx]
    splog, '12+log(O/H) Te vs KK04/R23: '+string(mindiv,format='(F4.2)')

    up = where((hii[indx].zt_log12oh_te gt mindiv),nup)
    if (nup ne 0L) then begin
       xup = hii[indx[up]].zt_log12oh_te
       xerrup = hii[indx[up]].zt_log12oh_te_err

       yup = hii[indx[up]].zstrong_12oh_kk04_r23_upper
       yerrup = hii[indx[up]].zstrong_12oh_kk04_r23_upper_err

       x = xup & xerr = xerrup
       y = yup & yerr = yerrup
    endif

    lo = where((hii[indx].zt_log12oh_te lt mindiv),nlo)
    if (nlo ne 0L) then begin
       xlo = hii[indx[lo]].zt_log12oh_te
       xerrlo = hii[indx[lo]].zt_log12oh_te_err

       ylo = hii[indx[lo]].zstrong_12oh_kk04_r23_lower
       yerrlo = hii[indx[lo]].zstrong_12oh_kk04_r23_lower_err

       x = [xup,xlo] & xerr = [xerrup,xerrlo]
       y = [yup,ylo] & yerr = [yerrup,yerrlo]
    endif

    upstats = im_stats(xup-yup)
    upxstr = strtrim(string(upstats.median,format='(F12.2)'),2)+'\pm'+$
      strtrim(string(upstats.sig68mean,format='(F12.2)'),2)

    lostats = im_stats(xlo-ylo)
    loxstr = strtrim(string(lostats.median,format='(F12.2)'),2)+'\pm'+$
      strtrim(string(lostats.sig68mean,format='(F12.2)'),2)

    stats = im_stats(x-y,/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)

    zindicators_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, talk=dotalk, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, charsize=charsize_6, $
      xstyle=3, ystyle=3, /right, /top, position=pos[*,2], hiipsize=0.8, /noerase
    djs_oplot, !x.crange, !y.crange, thick=postthick, color=djs_icolor(talkcolor), line=2
 
    legend, textoidl([upxstr+' (Upper)',loxstr+' (Lower)']), /right, /bottom, box=0, $
      charsize=charsize_5, charthick=postthick, textcolor=djs_icolor(talkcolor)
    legend, 'Kobulnicky & Kewley 2004', /left, /top, box=0, charsize=charsize_5, charthick=postthick, $
      textcolor=djs_icolor(talkcolor)

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; 12+log(O/H) Te vs 12+log(O/H) Empirical
; ------------------------------------------------------------

    psname = 'hii_12oh_te_vs_12oh_empirical_3panel'
    im_openclose, pspath+psname, postscript=postscript, encapsulated=encapsulated, xsize=8.5, ysize=10.5

    pagemaker, nx=1, ny=3, height=3.0*[1,1,1], width=6.8, xmargin=[1.3,0.4], $
      ymargin=[0.4,1.1], xspace=0, yspace=0, xpage=8.5, ypage=10.5, $
      position=pos, /normal

    if keyword_set(dotalk) then polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, $
      color=djs_icolor('black')

    ytitle = ohtitle
    xtitle = ohtitle+' T_{e}'

    xrange = ohrange
    yrange = ohrange
    
; --------------------------------------------------    
; [N II]/Ha
; --------------------------------------------------    

    indx = where((hii.zstrong_12oh_niiha_moustakas gt -900) and (hii.zt_log12oh_te gt -900.0),nindx)

    x = hii[indx].zt_log12oh_te
    xerr = hii[indx].zt_log12oh_te_err

    y = hii[indx].zstrong_12oh_niiha_moustakas
    yerr = hii[indx].zstrong_12oh_niiha_moustakas_err

    stats = im_stats(x-y,/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)
    
    zindicators_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, talk=dotalk, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=3, ystyle=3, /right, /top, position=pos[*,0], charsize=charsize_6, $
      xtickname=replicate(' ',10), hiipsize=0.8
    djs_oplot, !x.crange, !y.crange, thick=postthick, color=djs_icolor(talkcolor), line=2

    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_5, $
      charthick=postthick, textcolor=djs_icolor(talkcolor)
    legend, textoidl('[N II]/H\alpha'), /left, /top, box=0, $
      charsize=charsize_5, charthick=postthick, textcolor=djs_icolor(talkcolor)

; --------------------------------------------------    
; ([O III]/Hb)/([N II]/Ha)
; --------------------------------------------------    

    indx = where((hii.zstrong_12oh_oiiinii_moustakas gt -900) and (hii.zt_log12oh_te gt -900.0),nindx)

    x = hii[indx].zt_log12oh_te
    xerr = hii[indx].zt_log12oh_te_err

    y = hii[indx].zstrong_12oh_oiiinii_moustakas
    yerr = hii[indx].zstrong_12oh_oiiinii_moustakas_err

    stats = im_stats(x-y,/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)
    
    zindicators_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, talk=dotalk, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=3, ystyle=3, /right, /top, position=pos[*,1], charsize=charsize_6, $
      xtickname=replicate(' ',10), hiipsize=0.8, /noerase
    djs_oplot, !x.crange, !y.crange, thick=postthick, color=djs_icolor(talkcolor), line=2

    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_5, $
      charthick=postthick, textcolor=djs_icolor(talkcolor)
    legend, textoidl('([O III]/H\beta)/([N II]/H\alpha)'), /left, /top, box=0, $
      charsize=charsize_5, charthick=postthick, textcolor=djs_icolor(talkcolor)

; --------------------------------------------------    
; P-method
; --------------------------------------------------    

    indx = where((hii.zt_log12oh_te gt -900.0) and (hii.zstrong_12oh_p01_upper gt -900.0) and $
      (hii.zstrong_12oh_p01_lower gt -900.0),nindx)

    up = where((hii[indx].zt_log12oh_te gt 8.2) and (hii[indx].zstrong_12oh_p01_upper gt 8.2),nup)
    if (nup ne 0L) then begin
       xup = hii[indx[up]].zt_log12oh_te
       xerrup = hii[indx[up]].zt_log12oh_te_err

       yup = hii[indx[up]].zstrong_12oh_p01_upper
       yerrup = hii[indx[up]].zstrong_12oh_p01_upper_err

       x = xup & xerr = xerrup
       y = yup & yerr = yerrup
    endif

    lo = where((hii[indx].zt_log12oh_te lt 8.0),nlo)
;   lo = where((hii[indx].zt_log12oh_te lt 8.0) and (hii[indx].zstrong_12oh_p01_lower lt 8.0),nlo)
    if (nlo ne 0L) then begin
       xlo = hii[indx[lo]].zt_log12oh_te
       xerrlo = hii[indx[lo]].zt_log12oh_te_err

       ylo = hii[indx[lo]].zstrong_12oh_p01_lower
       yerrlo = hii[indx[lo]].zstrong_12oh_p01_lower_err

       x = [xup,xlo] & xerr = [xerrup,xerrlo]
       y = [yup,ylo] & yerr = [yerrup,yerrlo]
    endif

    stats = im_stats(x-y,/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)

    zindicators_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, talk=dotalk, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=3, ystyle=3, /right, /top, position=pos[*,2], charsize=charsize_6, $
      hiipsize=0.8, /noerase
    djs_oplot, !x.crange, !y.crange, thick=postthick, color=djs_icolor(talkcolor), line=2

    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_5, $
      charthick=postthick, textcolor=djs_icolor(talkcolor)
    legend, textoidl('P-method'), /left, /top, box=0, $
      charsize=charsize_5, charthick=postthick, textcolor=djs_icolor(talkcolor)

    im_openclose, postscript=postscript, /close    
    
; ###########################################################################
; line ratios versus abundance
; ###########################################################################

; HERE!    
    
; ###########################################################################
; abundance versus line ratios
; ###########################################################################

; ------------------------------------------------------------
; 12+log (O/H) [Te] vs [O II]/Ha
; ------------------------------------------------------------

    psname = 'hii_12oh_te_vs_oiiha'
    im_openclose, pspath+psname, postscript=postscript, encapsulated=encapsulated, xsize=8.5, ysize=8.0

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.6, height=6.6, $
      xmargin=[1.5,0.4], ymargin=[0.3,1.1], xpage=8.5, ypage=8.0, $
      position=pos, /normal

    if keyword_set(dotalk) then polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, $
      color=djs_icolor('black')
    
    indx = where((hii.oii_h_alpha gt -900.0) and (hii.zt_log12oh_te gt -900.0),nindx)
    
    x = hii[indx].zt_log12oh_te
    xerr = hii[indx].zt_log12oh_te_err

    y = hii[indx].oii_h_alpha
    yerr = hii[indx].oii_h_alpha_err

;   xrange = [7.0,9.4]
    xrange = ohrange
    yrange = oiiharange

    xtitle = '12 + log (O/H) [T_{e}]'
    ytitle = 'log ([O II] \lambda3727 / H\alpha)'

    zindicators_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, talk=dotalk, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=3, ystyle=3, /right, /top, hiipsize=1.0, position=pos[*,0]
    
; overlay some representative kewley grids

    hiimodel1 = interpol(hiigrids[*,1].oii_h_alpha,log12ohgrid,log12ohplot,/quad)
    hiimodel2 = interpol(hiigrids[*,5].oii_h_alpha,log12ohgrid,log12ohplot,/quad)
    sbmodel1 = interpol(sbgrids[*,1].oii_h_alpha,log12ohgrid,log12ohplot,/quad)
    sbmodel2 = interpol(sbgrids[*,5].oii_h_alpha,log12ohgrid,log12ohplot,/quad)
    
    djs_oplot, log12ohplot, hiimodel1, line=2, thick=postthick, color=djs_icolor(talkcolor)
    djs_oplot, log12ohplot, hiimodel2, line=2, thick=postthick, color=djs_icolor(talkcolor)
    djs_oplot, log12ohplot, sbmodel1, line=0, thick=postthick, color=djs_icolor(talkcolor)
    djs_oplot, log12ohplot, sbmodel2, line=0, thick=postthick, color=djs_icolor(talkcolor)
    
    xyouts, log12ohplot[0], hiimodel1[0], string(U[1],format='(F5.2)'), /data, charsize=1.5, $
      charthick=postthick, align=1.0, color=djs_icolor(talkcolor)
    xyouts, log12ohplot[0], hiimodel2[0], string(U[5],format='(F5.2)'), /data, charsize=1.5, $
      charthick=postthick, align=1.0, color=djs_icolor(talkcolor)

    legend, ['HII Region','Starburst'], line=[2,0], /left, /top, box=0, color=djs_icolor(talkcolor), $
      charsize=1.5, charthick=postthick, thick=postthick, textcolor=djs_icolor(talkcolor)

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; 12+log (O/H) [Te] vs [N II]/Ha
; ------------------------------------------------------------

    psname = 'hii_12oh_te_vs_niiha'
    im_openclose, pspath+psname, postscript=postscript, encapsulated=encapsulated, xsize=8.5, ysize=8.0

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.6, height=6.6, $
      xmargin=[1.5,0.4], ymargin=[0.3,1.1], xpage=8.5, ypage=8.0, $
      position=pos, /normal

    if keyword_set(dotalk) then polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, $
      color=djs_icolor('black')

    indx = where((hii.nii_6584_h_alpha gt -900.0) and (hii.zt_log12oh_te gt -900.0),nindx)
    
    x = hii[indx].zt_log12oh_te
    xerr = hii[indx].zt_log12oh_te_err

    y = hii[indx].nii_6584_h_alpha
    yerr = hii[indx].nii_6584_h_alpha_err

    xrange = ohrange
    yrange = niiharange

    xtitle = '12 + log (O/H) [T_{e}]'
    ytitle = 'log ([N II] \lambda6584 / H\alpha)'
    
    zindicators_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, talk=dotalk, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /left, xstyle=3, ystyle=3, /top, hiipsize=1.0, position=pos[*,0]

; overlay some representative kewley grids

    hiimodel1 = interpol(hiigrids[*,1].nii_6584_h_alpha,log12ohgrid,log12ohplot,/quad)
    hiimodel2 = interpol(hiigrids[*,5].nii_6584_h_alpha,log12ohgrid,log12ohplot,/quad)
    sbmodel1 = interpol(sbgrids[*,1].nii_6584_h_alpha,log12ohgrid,log12ohplot,/quad)
    sbmodel2 = interpol(sbgrids[*,5].nii_6584_h_alpha,log12ohgrid,log12ohplot,/quad)
    
    djs_oplot, log12ohplot, hiimodel1, line=2, thick=postthick, color=djs_icolor(talkcolor)
    djs_oplot, log12ohplot, hiimodel2, line=2, thick=postthick, color=djs_icolor(talkcolor)
    djs_oplot, log12ohplot, sbmodel1, line=0, thick=postthick, color=djs_icolor(talkcolor)
    djs_oplot, log12ohplot, sbmodel2, line=0, thick=postthick, color=djs_icolor(talkcolor)
    
    xyouts, log12ohplot[0], hiimodel1[0], string(U[1],format='(F5.2)'), /data, charsize=1.5, $
      charthick=postthick, align=1.0, color=djs_icolor(talkcolor)
    xyouts, log12ohplot[0], hiimodel2[0], string(U[5],format='(F5.2)'), /data, charsize=1.5, $
      charthick=postthick, align=1.0, color=djs_icolor(talkcolor)
    
    legend, ['HII Region','Starburst'], line=[2,0], /left, /top, box=0, color=djs_icolor(talkcolor), $
      charsize=1.5, charthick=postthick, thick=postthick, textcolor=djs_icolor(talkcolor)

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; 12+log (O/H) [Te] vs [N II]/[O II]
; ------------------------------------------------------------

    psname = 'hii_12oh_te_vs_niioii'
    im_openclose, pspath+psname, postscript=postscript, encapsulated=encapsulated, xsize=8.5, ysize=8.0

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.6, height=6.6, $
      xmargin=[1.5,0.4], ymargin=[0.3,1.1], xpage=8.5, ypage=8.0, $
      position=pos, /normal

    if keyword_set(dotalk) then polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, $
      color=djs_icolor('black')

    indx = where((hii.nii_6584_oii gt -900.0) and (hii.zt_log12oh_te gt -900.0),nindx)
    
    x = hii[indx].zt_log12oh_te
    xerr = hii[indx].zt_log12oh_te_err

    y = hii[indx].nii_6584_oii
    yerr = hii[indx].nii_6584_oii_err

    xrange = ohrange
    yrange = niioiirange

    xtitle = '12 + log (O/H) [T_{e}]'
    ytitle = 'log ([N II] \lambda6584 / [O II] \lambda3727)'
    
    zindicators_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, talk=dotalk, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=3, ystyle=3, hiipsize=1.0, position=pos[*,0]

; overlay the kewley grids

    hiimodel1 = interpol(hiigrids[*,1].nii_6584_oii,log12ohgrid,log12ohplot,/quad)
    hiimodel2 = interpol(hiigrids[*,5].nii_6584_oii,log12ohgrid,log12ohplot,/quad)
    sbmodel1 = interpol(sbgrids[*,1].nii_6584_oii,log12ohgrid,log12ohplot,/quad)
    sbmodel2 = interpol(sbgrids[*,5].nii_6584_oii,log12ohgrid,log12ohplot,/quad)
    
    djs_oplot, log12ohplot, hiimodel1, line=2, thick=postthick, color=djs_icolor(talkcolor)
    djs_oplot, log12ohplot, hiimodel2, line=2, thick=postthick, color=djs_icolor(talkcolor)
    djs_oplot, log12ohplot, sbmodel1, line=0, thick=postthick, color=djs_icolor(talkcolor)
    djs_oplot, log12ohplot, sbmodel2, line=0, thick=postthick, color=djs_icolor(talkcolor)
    
    xyouts, log12ohplot[0], hiimodel1[0], string(U[1],format='(F5.2)'), /data, charsize=1.5, $
      charthick=postthick, align=1.0, color=djs_icolor(talkcolor)
    xyouts, log12ohplot[0], hiimodel2[0], string(U[5],format='(F5.2)'), /data, charsize=1.5, $
      charthick=postthick, align=1.0, color=djs_icolor(talkcolor)

    legend, ['HII Region','Starburst'], line=[2,0], /left, /top, box=0, color=djs_icolor(talkcolor), $
      charsize=1.5, charthick=postthick, thick=postthick, textcolor=djs_icolor(talkcolor)
       
    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; 12+log (O/H) [Te] vs [O III]/Hb
; ------------------------------------------------------------

    psname = 'hii_12oh_te_vs_oiiihb'
    im_openclose, pspath+psname, postscript=postscript, encapsulated=encapsulated, xsize=8.5, ysize=8.0

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.6, height=6.6, $
      xmargin=[1.5,0.4], ymargin=[0.3,1.1], xpage=8.5, ypage=8.0, $
      position=pos, /normal

    if keyword_set(dotalk) then polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, $
      color=djs_icolor('black')

    indx = where((hii.oiii_5007_h_beta gt -900.0) and (hii.zt_log12oh_te gt -900.0),nindx)
    
    x = hii[indx].zt_log12oh_te
    xerr = hii[indx].zt_log12oh_te_err

    y = hii[indx].oiii_5007_h_beta
    yerr = hii[indx].oiii_5007_h_beta_err

    xrange = ohrange
    yrange = oiiihbrange

    xtitle = '12 + log (O/H) [T_{e}]'
    ytitle = 'log ([O III] \lambda5007 / H\beta)'

    zindicators_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, talk=dotalk, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=3, ystyle=3, /right, /top, hiipsize=1.0, position=pos[*,0]

; overlay some representative kewley grids

    hiimodel1 = interpol(hiigrids[*,1].oiii_5007_h_beta,log12ohgrid,log12ohplot,/quad)
    hiimodel2 = interpol(hiigrids[*,5].oiii_5007_h_beta,log12ohgrid,log12ohplot,/quad)
    sbmodel1 = interpol(sbgrids[*,1].oiii_5007_h_beta,log12ohgrid,log12ohplot,/quad)
    sbmodel2 = interpol(sbgrids[*,5].oiii_5007_h_beta,log12ohgrid,log12ohplot,/quad)
    
    djs_oplot, log12ohplot, hiimodel1, line=2, thick=postthick, color=djs_icolor(talkcolor)
    djs_oplot, log12ohplot, hiimodel2, line=2, thick=postthick, color=djs_icolor(talkcolor)
    djs_oplot, log12ohplot, sbmodel1, line=0, thick=postthick, color=djs_icolor(talkcolor)
    djs_oplot, log12ohplot, sbmodel2, line=0, thick=postthick, color=djs_icolor(talkcolor)
    
    xyouts, log12ohplot[0], hiimodel1[0], string(U[1],format='(F5.2)'), /data, charsize=1.5, $
      charthick=postthick, align=1.0, color=djs_icolor(talkcolor)
    xyouts, log12ohplot[0], hiimodel2[0], string(U[5],format='(F5.2)'), /data, charsize=1.5, $
      charthick=postthick, align=1.0, color=djs_icolor(talkcolor)

    legend, ['HII Region','Starburst'], line=[2,0], /left, /top, box=0, color=djs_icolor(talkcolor), $
      charsize=1.5, charthick=postthick, thick=postthick, textcolor=djs_icolor(talkcolor)

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; 12+log (O/H) [Te] vs [N II]/[S II]
; ------------------------------------------------------------

    psname = 'hii_12oh_te_vs_niisii'
    im_openclose, pspath+psname, postscript=postscript, encapsulated=encapsulated, xsize=8.5, ysize=8.0

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.6, height=6.6, $
      xmargin=[1.5,0.4], ymargin=[0.3,1.1], xpage=8.5, ypage=8.0, $
      position=pos, /normal

    if keyword_set(dotalk) then polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, $
      color=djs_icolor('black')

    indx = where((hii.nii_6584_sii gt -900.0) and (hii.zt_log12oh_te gt -900.0),nindx)
    
    x = hii[indx].zt_log12oh_te
    xerr = hii[indx].zt_log12oh_te_err

    y = hii[indx].nii_6584_sii
    yerr = hii[indx].nii_6584_sii_err

    xrange = ohrange
    yrange = niisiirange

    xtitle = '12 + log (O/H) [T_{e}]'
    ytitle = 'log ([N II] \lambda6584 / [S II] \lambda\lambda6716,6731)'

    zindicators_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, talk=dotalk, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /left, xstyle=3, ystyle=3, /right, /top, hiipsize=1.0, position=pos[*,0]

; overlay some representative kewley grids

    hiimodel1 = interpol(hiigrids[*,1].nii_6584_sii,log12ohgrid,log12ohplot,/quad)
    hiimodel2 = interpol(hiigrids[*,5].nii_6584_sii,log12ohgrid,log12ohplot,/quad)
    sbmodel1 = interpol(sbgrids[*,1].nii_6584_sii,log12ohgrid,log12ohplot,/quad)
    sbmodel2 = interpol(sbgrids[*,5].nii_6584_sii,log12ohgrid,log12ohplot,/quad)
    
    djs_oplot, log12ohplot, hiimodel1, line=2, thick=postthick, color=djs_icolor(talkcolor)
    djs_oplot, log12ohplot, hiimodel2, line=2, thick=postthick, color=djs_icolor(talkcolor)
    djs_oplot, log12ohplot, sbmodel1, line=0, thick=postthick, color=djs_icolor(talkcolor)
    djs_oplot, log12ohplot, sbmodel2, line=0, thick=postthick, color=djs_icolor(talkcolor)
    
    xyouts, log12ohplot[0], hiimodel1[0], string(U[1],format='(F5.2)'), /data, charsize=1.5, $
      charthick=postthick, align=1.0, color=djs_icolor(talkcolor)
    xyouts, log12ohplot[0], hiimodel2[0], string(U[5],format='(F5.2)'), /data, charsize=1.5, $
      charthick=postthick, align=1.0, color=djs_icolor(talkcolor)

    legend, ['HII Region','Starburst'], line=[2,0], /left, /top, box=0, color=djs_icolor(talkcolor), $
      charsize=1.5, charthick=postthick, thick=postthick, textcolor=djs_icolor(talkcolor)
       
    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; 12+log (O/H) [Te] vs [O II]/[S II]
; ------------------------------------------------------------

    psname = 'hii_12oh_te_vs_oiisii'
    im_openclose, pspath+psname, postscript=postscript, encapsulated=encapsulated, xsize=8.5, ysize=8.0

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.6, height=6.6, $
      xmargin=[1.5,0.4], ymargin=[0.3,1.1], xpage=8.5, ypage=8.0, $
      position=pos, /normal

    if keyword_set(dotalk) then polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, $
      color=djs_icolor('black')

    indx = where((hii.oii_sii gt -900.0) and (hii.zt_log12oh_te gt -900.0),nindx)
    
    x = hii[indx].zt_log12oh_te
    xerr = hii[indx].zt_log12oh_te_err

    y = hii[indx].oii_sii
    yerr = hii[indx].oii_sii_err

    xrange = ohrange
    yrange = oiisiirange

    xtitle = '12 + log (O/H) [T_{e}]'
    ytitle = 'log ([O II] \lambda3727 / [S II] \lambda\lambda6716,6731)'

    zindicators_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, talk=dotalk, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /left, xstyle=3, ystyle=3, /right, /top, hiipsize=1.0, position=pos[*,0]

; overlay some representative kewley grids

    hiimodel1 = interpol(hiigrids[*,1].oii_sii,log12ohgrid,log12ohplot,/quad)
    hiimodel2 = interpol(hiigrids[*,5].oii_sii,log12ohgrid,log12ohplot,/quad)
    sbmodel1 = interpol(sbgrids[*,1].oii_sii,log12ohgrid,log12ohplot,/quad)
    sbmodel2 = interpol(sbgrids[*,5].oii_sii,log12ohgrid,log12ohplot,/quad)
    
    djs_oplot, log12ohplot, hiimodel1, line=2, thick=postthick, color=djs_icolor(talkcolor)
    djs_oplot, log12ohplot, hiimodel2, line=2, thick=postthick, color=djs_icolor(talkcolor)
    djs_oplot, log12ohplot, sbmodel1, line=0, thick=postthick, color=djs_icolor(talkcolor)
    djs_oplot, log12ohplot, sbmodel2, line=0, thick=postthick, color=djs_icolor(talkcolor)
    
    xyouts, log12ohplot[0], hiimodel1[0], string(U[1],format='(F5.2)'), /data, charsize=1.5, $
      charthick=postthick, align=1.0, color=djs_icolor(talkcolor)
    xyouts, log12ohplot[0], hiimodel2[0], string(U[5],format='(F5.2)'), /data, charsize=1.5, $
      charthick=postthick, align=1.0, color=djs_icolor(talkcolor)

    legend, ['HII Region','Starburst'], line=[2,0], /left, /bottom, box=0, color=djs_icolor(talkcolor), $
      charsize=1.5, charthick=postthick, thick=postthick, textcolor=djs_icolor(talkcolor)
       
    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; 12+log (O/H) [Te] vs [N II]/[O III]
; ------------------------------------------------------------

    psname = 'hii_12oh_te_vs_niioiii'
    im_openclose, pspath+psname, postscript=postscript, encapsulated=encapsulated, xsize=8.5, ysize=8.0

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.6, height=6.6, $
      xmargin=[1.5,0.4], ymargin=[0.3,1.1], xpage=8.5, ypage=8.0, $
      position=pos, /normal

    if keyword_set(dotalk) then polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, $
      color=djs_icolor('black')

    indx = where((hii.nii_6584_oiii_5007 gt -900.0) and (hii.zt_log12oh_te gt -900.0),nindx)
    
    x = hii[indx].zt_log12oh_te
    xerr = hii[indx].zt_log12oh_te_err

    y = hii[indx].nii_6584_oiii_5007
    yerr = hii[indx].nii_6584_oiii_5007_err

    xrange = ohrange
    yrange = niioiiirange
    
    xtitle = '12 + log (O/H) [T_{e}]'
    ytitle = 'log ([N II] \lambda6584 / [O III] \lambda5007)'
    
    zindicators_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, talk=dotalk, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /left, xstyle=3, ystyle=3, /right, /top, hiipsize=1.0, position=pos[*,0]

; overlay some representative kewley grids

    hiimodel1 = interpol(hiigrids[*,1].nii_6584_oiii_5007,log12ohgrid,log12ohplot,/quad)
    hiimodel2 = interpol(hiigrids[*,5].nii_6584_oiii_5007,log12ohgrid,log12ohplot,/quad)
    sbmodel1 = interpol(sbgrids[*,1].nii_6584_oiii_5007,log12ohgrid,log12ohplot,/quad)
    sbmodel2 = interpol(sbgrids[*,5].nii_6584_oiii_5007,log12ohgrid,log12ohplot,/quad)
    
    djs_oplot, log12ohplot, hiimodel1, line=2, thick=postthick, color=djs_icolor(talkcolor)
    djs_oplot, log12ohplot, hiimodel2, line=2, thick=postthick, color=djs_icolor(talkcolor)
    djs_oplot, log12ohplot, sbmodel1, line=0, thick=postthick, color=djs_icolor(talkcolor)
    djs_oplot, log12ohplot, sbmodel2, line=0, thick=postthick, color=djs_icolor(talkcolor)
    
    xyouts, log12ohplot[0], hiimodel1[0], string(U[1],format='(F5.2)'), /data, charsize=1.5, $
      charthick=postthick, align=1.0, color=djs_icolor(talkcolor)
    xyouts, log12ohplot[0], hiimodel2[0], string(U[5],format='(F5.2)'), /data, charsize=1.5, $
      charthick=postthick, align=1.0, color=djs_icolor(talkcolor)

    legend, ['HII Region','Starburst'], line=[2,0], /left, /top, box=0, color=djs_icolor(talkcolor), $
      charsize=1.5, charthick=postthick, thick=postthick, textcolor=djs_icolor(talkcolor)
       
    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; 12+log (O/H) [Te] vs [O III]/[O II]
; ------------------------------------------------------------

    psname = 'hii_12oh_te_vs_oiiioii'
    im_openclose, pspath+psname, postscript=postscript, encapsulated=encapsulated, xsize=8.5, ysize=8.0

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.6, height=6.6, $
      xmargin=[1.5,0.4], ymargin=[0.3,1.1], xpage=8.5, ypage=8.0, $
      position=pos, /normal

    if keyword_set(dotalk) then polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, $
      color=djs_icolor('black')

    indx = where((hii.oiii_5007_oii gt -900.0) and (hii.zt_log12oh_te gt -900.0),nindx)
    
    x = hii[indx].zt_log12oh_te
    xerr = hii[indx].zt_log12oh_te_err

    y = hii[indx].oiii_5007_oii
    yerr = hii[indx].oiii_5007_oii_err

    xrange = ohrange
    yrange = oiiioiirange

    xtitle = '12 + log (O/H) [T_{e}]'
    ytitle = 'log ([O III] \lambda5007 / [O II] \lambda3727)'
    
    zindicators_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, talk=dotalk, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /left, xstyle=3, ystyle=3, /right, /top, hiipsize=1.0, position=pos[*,0]

; overlay some representative kewley grids

    hiimodel1 = interpol(hiigrids[*,1].oiii_5007_oii,log12ohgrid,log12ohplot,/quad)
    hiimodel2 = interpol(hiigrids[*,5].oiii_5007_oii,log12ohgrid,log12ohplot,/quad)
    sbmodel1 = interpol(sbgrids[*,1].oiii_5007_oii,log12ohgrid,log12ohplot,/quad)
    sbmodel2 = interpol(sbgrids[*,5].oiii_5007_oii,log12ohgrid,log12ohplot,/quad)
    
    djs_oplot, log12ohplot, hiimodel1, line=2, thick=postthick, color=djs_icolor(talkcolor)
    djs_oplot, log12ohplot, hiimodel2, line=2, thick=postthick, color=djs_icolor(talkcolor)
    djs_oplot, log12ohplot, sbmodel1, line=0, thick=postthick, color=djs_icolor(talkcolor)
    djs_oplot, log12ohplot, sbmodel2, line=0, thick=postthick, color=djs_icolor(talkcolor)
    
    xyouts, log12ohplot[0], hiimodel1[0], string(U[1],format='(F5.2)'), /data, $
      charsize=1.5, charthick=postthick, align=1.0, color=djs_icolor(talkcolor)
    xyouts, log12ohplot[0], hiimodel2[0], string(U[5],format='(F5.2)'), /data, $
      charsize=1.5, charthick=postthick, align=1.0, color=djs_icolor(talkcolor)

    legend, ['HII Region','Starburst'], line=[2,0], /right, /top, box=0, color=djs_icolor(talkcolor), $
      charsize=1.5, charthick=postthick, thick=postthick, textcolor=djs_icolor(talkcolor)
       
    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; 12+log (O/H) [Te] vs P
; ------------------------------------------------------------

    psname = 'hii_12oh_te_vs_P'
    im_openclose, pspath+psname, postscript=postscript, encapsulated=encapsulated, xsize=8.5, ysize=8.0

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.6, height=6.6, $
      xmargin=[1.5,0.4], ymargin=[0.3,1.1], xpage=8.5, ypage=8.0, $
      position=pos, /normal

    if keyword_set(dotalk) then polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, $
      color=djs_icolor('black')

    indx = where((hii.zstrong_P gt -900.0) and (hii.zt_log12oh_te gt -900.0),nindx)
    
    x = hii[indx].zt_log12oh_te
    xerr = hii[indx].zt_log12oh_te_err

    y = hii[indx].zstrong_P
    yerr = hii[indx].zstrong_P_err

    xrange = ohrange
    yrange = prange

    xtitle = '12 + log (O/H) [T_{e}]'
    ytitle = 'P'
    
    zindicators_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, talk=dotalk, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /left, xstyle=3, ystyle=3, /right, /top, hiipsize=1.0, position=pos[*,0]

; overlay some representative kewley grids

    hiimodel1 = interpol(hiigrids[*,1].P,log12ohgrid,log12ohplot,/quad)
    hiimodel2 = interpol(hiigrids[*,5].P,log12ohgrid,log12ohplot,/quad)
    sbmodel1 = interpol(sbgrids[*,1].P,log12ohgrid,log12ohplot,/quad)
    sbmodel2 = interpol(sbgrids[*,5].P,log12ohgrid,log12ohplot,/quad)
    
    djs_oplot, log12ohplot, hiimodel1, line=2, thick=postthick, color=djs_icolor(talkcolor)
    djs_oplot, log12ohplot, hiimodel2, line=2, thick=postthick, color=djs_icolor(talkcolor)
    djs_oplot, log12ohplot, sbmodel1, line=0, thick=postthick, color=djs_icolor(talkcolor)
    djs_oplot, log12ohplot, sbmodel2, line=0, thick=postthick, color=djs_icolor(talkcolor)
    
    xyouts, log12ohplot[0], hiimodel1[0], string(U[1],format='(F5.2)'), /data, $
      charsize=1.5, charthick=postthick, align=1.0, color=djs_icolor(talkcolor)
    xyouts, log12ohplot[0], hiimodel2[0], string(U[5],format='(F5.2)'), /data, $
      charsize=1.5, charthick=postthick, align=1.0, color=djs_icolor(talkcolor)

    legend, ['HII Region','Starburst'], line=[2,0], /right, /top, box=0, color=djs_icolor(talkcolor), $
      charsize=1.5, charthick=postthick, thick=postthick, textcolor=djs_icolor(talkcolor)
       
    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; 12+log (O/H) [Te] vs [S II]/Ha
; ------------------------------------------------------------

    psname = 'hii_12oh_te_vs_siiha'
    im_openclose, pspath+psname, postscript=postscript, encapsulated=encapsulated, xsize=8.5, ysize=8.0

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.6, height=6.6, $
      xmargin=[1.5,0.4], ymargin=[0.3,1.1], xpage=8.5, ypage=8.0, $
      position=pos, /normal

    if keyword_set(dotalk) then polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, $
      color=djs_icolor('black')

    indx = where((hii.sii_h_alpha gt -900.0) and (hii.zt_log12oh_te gt -900.0),nindx)
    
    x = hii[indx].zt_log12oh_te
    xerr = hii[indx].zt_log12oh_te_err

    y = hii[indx].sii_h_alpha
    yerr = hii[indx].sii_h_alpha_err

    xrange = ohrange
    yrange = siiharange
    
    xtitle = '12 + log (O/H) [T_{e}]'
    ytitle = 'log ([S II] \lambda\lambda6716,6731 / H\alpha)'
    
    zindicators_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, talk=dotalk, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /left, xstyle=3, ystyle=3, /right, /top, hiipsize=1.0, position=pos[*,0]

; overlay some representative kewley grids

    hiimodel1 = interpol(hiigrids[*,1].sii_h_alpha,log12ohgrid,log12ohplot,/quad)
    hiimodel2 = interpol(hiigrids[*,5].sii_h_alpha,log12ohgrid,log12ohplot,/quad)
    sbmodel1 = interpol(sbgrids[*,1].sii_h_alpha,log12ohgrid,log12ohplot,/quad)
    sbmodel2 = interpol(sbgrids[*,5].sii_h_alpha,log12ohgrid,log12ohplot,/quad)
    
    djs_oplot, log12ohplot, hiimodel1, line=2, thick=postthick, color=djs_icolor(talkcolor)
    djs_oplot, log12ohplot, hiimodel2, line=2, thick=postthick, color=djs_icolor(talkcolor)
    djs_oplot, log12ohplot, sbmodel1, line=0, thick=postthick, color=djs_icolor(talkcolor)
    djs_oplot, log12ohplot, sbmodel2, line=0, thick=postthick, color=djs_icolor(talkcolor)
    
    xyouts, log12ohplot[0], hiimodel1[0], string(U[1],format='(F5.2)'), /data, charsize=1.5, $
      charthick=postthick, align=1.0, color=djs_icolor(talkcolor)
    xyouts, log12ohplot[0], hiimodel2[0], string(U[5],format='(F5.2)'), /data, charsize=1.5, $
      charthick=postthick, align=1.0, color=djs_icolor(talkcolor)

    legend, ['HII Region','Starburst'], line=[2,0], /left, /top, box=0, color=djs_icolor(talkcolor), $
      charsize=1.5, charthick=postthick, thick=postthick, textcolor=djs_icolor(talkcolor)
       
    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; 12+log (O/H) [Te] vs ([O III]/Hb)/([N II]/Ha) - Pettini & Pagel (2004)
; ------------------------------------------------------------

    psname = 'hii_12oh_te_vs_oiiinii'
    im_openclose, pspath+psname, postscript=postscript, encapsulated=encapsulated, xsize=8.5, ysize=8.0

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.6, height=6.6, $
      xmargin=[1.5,0.4], ymargin=[0.3,1.1], xpage=8.5, ypage=8.0, $
      position=pos, /normal

    if keyword_set(dotalk) then polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, $
      color=djs_icolor('black')

    indx = where((hii.oiii_5007_h_beta_nii_6584_h_alpha gt -900.0) and (hii.zt_log12oh_te gt -900.0),nindx)
    
    x = hii[indx].zt_log12oh_te
    xerr = hii[indx].zt_log12oh_te_err

    y = hii[indx].oiii_5007_h_beta_nii_6584_h_alpha
    yerr = hii[indx].oiii_5007_h_beta_nii_6584_h_alpha_err

    xrange = ohrange
    yrange = o3n2range
    
    xtitle = '12 + log (O/H) [T_{e}]'
    ytitle = 'log {([O III] \lambda5007/H\beta)/([N II] \lambda6584/H\alpha)}'
    
    zindicators_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, talk=dotalk, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /left, xstyle=3, ystyle=3, /right, /top, hiipsize=1.0, position=pos[*,0]

; overlay some representative kewley grids

    hiimodel1 = interpol(hiigrids[*,1].oiii_5007_h_beta_nii_6584_h_alpha,log12ohgrid,log12ohplot,/quad)
    hiimodel2 = interpol(hiigrids[*,5].oiii_5007_h_beta_nii_6584_h_alpha,log12ohgrid,log12ohplot,/quad)
    sbmodel1 = interpol(sbgrids[*,1].oiii_5007_h_beta_nii_6584_h_alpha,log12ohgrid,log12ohplot,/quad)
    sbmodel2 = interpol(sbgrids[*,5].oiii_5007_h_beta_nii_6584_h_alpha,log12ohgrid,log12ohplot,/quad)
    
    djs_oplot, log12ohplot, hiimodel1, line=2, thick=postthick, color=djs_icolor(talkcolor)
    djs_oplot, log12ohplot, hiimodel2, line=2, thick=postthick, color=djs_icolor(talkcolor)
    djs_oplot, log12ohplot, sbmodel1, line=0, thick=postthick, color=djs_icolor(talkcolor)
    djs_oplot, log12ohplot, sbmodel2, line=0, thick=postthick, color=djs_icolor(talkcolor)
    
    xyouts, log12ohplot[0], hiimodel1[0], string(U[1],format='(F5.2)'), /data, charsize=1.5, $
      charthick=postthick, align=1.0, color=djs_icolor(talkcolor)
    xyouts, log12ohplot[0], hiimodel2[0], string(U[5],format='(F5.2)'), /data, charsize=1.5, $
      charthick=postthick, align=1.0, color=djs_icolor(talkcolor)

    legend, ['HII Region','Starburst'], line=[2,0], /right, /top, box=0, color=djs_icolor(talkcolor), $
      charsize=1.5, charthick=postthick, thick=postthick, textcolor=djs_icolor(talkcolor)

    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; 12+log (O/H) [Te] vs ([O II]/Hb)*([O II]/[S II]): Dopita & Evans
; (1986) 
; ------------------------------------------------------------

    psname = 'hii_12oh_te_vs_oiihb_oiisii'
    im_openclose, pspath+psname, postscript=postscript, encapsulated=encapsulated, xsize=8.5, ysize=8.0

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.6, height=6.6, $
      xmargin=[1.5,0.4], ymargin=[0.3,1.1], xpage=8.5, ypage=8.0, $
      position=pos, /normal

    if keyword_set(dotalk) then polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, $
      color=djs_icolor('black')

    indx = where((hii.oii_h_beta_oii_sii gt -900.0) and (hii.zt_log12oh_te gt -900.0),nindx)
    
    x = hii[indx].zt_log12oh_te
    xerr = hii[indx].zt_log12oh_te_err

    y = hii[indx].oii_h_beta_oii_sii
    yerr = hii[indx].oii_h_beta_oii_sii_err

    xrange = ohrange
    yrange = oiihboiisiirange

    xtitle = '12 + log (O/H) [T_{e}]'
    ytitle = 'log {([O II] \lambda3727/H\beta) ([O II] \lambda3727/ [S II] \lambda\lambda6716,6731)}'
    
    zindicators_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, talk=dotalk, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /left, xstyle=3, ystyle=3, /right, /top, charsize=1.5, hiipsize=1.0, position=pos[*,0]

; overlay some representative kewley grids

    hiimodel1 = interpol(hiigrids[*,1].oii_h_beta_oii_sii,log12ohgrid,log12ohplot,/quad)
    hiimodel2 = interpol(hiigrids[*,5].oii_h_beta_oii_sii,log12ohgrid,log12ohplot,/quad)
    sbmodel1 = interpol(sbgrids[*,1].oii_h_beta_oii_sii,log12ohgrid,log12ohplot,/quad)
    sbmodel2 = interpol(sbgrids[*,5].oii_h_beta_oii_sii,log12ohgrid,log12ohplot,/quad)
    
    djs_oplot, log12ohplot, hiimodel1, line=2, thick=postthick, color=djs_icolor(talkcolor)
    djs_oplot, log12ohplot, hiimodel2, line=2, thick=postthick, color=djs_icolor(talkcolor)
    djs_oplot, log12ohplot, sbmodel1, line=0, thick=postthick, color=djs_icolor(talkcolor)
    djs_oplot, log12ohplot, sbmodel2, line=0, thick=postthick, color=djs_icolor(talkcolor)
    
    xyouts, log12ohplot[0], hiimodel1[0], string(U[1],format='(F5.2)'), /data, charsize=1.5, $
      charthick=postthick, align=1.0, color=djs_icolor(talkcolor)
    xyouts, log12ohplot[0], hiimodel2[0], string(U[5],format='(F5.2)'), /data, charsize=1.5, $
      charthick=postthick, align=1.0, color=djs_icolor(talkcolor)

    legend, ['HII Region','Starburst'], line=[2,0], /right, /top, box=0, color=djs_icolor(talkcolor), $
      charsize=1.5, charthick=postthick, thick=postthick, textcolor=djs_icolor(talkcolor)
       
    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; 12+log (O/H) [Te] vs r23
; ------------------------------------------------------------

    psname = 'hii_12oh_te_vs_r23'
    im_openclose, pspath+psname, postscript=postscript, encapsulated=encapsulated, xsize=8.5, ysize=8.0

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.6, height=6.6, $
      xmargin=[1.5,0.4], ymargin=[0.3,1.1], xpage=8.5, ypage=8.0, $
      position=pos, /normal

    if keyword_set(dotalk) then polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, $
      color=djs_icolor('black')

    indx = where((hii.zstrong_r23 gt -900.0) and (hii.zt_log12oh_te gt -900.0),nindx)
    
    x = hii[indx].zt_log12oh_te
    xerr = hii[indx].zt_log12oh_te_err

    y = alog10(hii[indx].zstrong_r23)
    yerr = hii[indx].zstrong_r23_err/hii[indx].zstrong_r23/alog(10.0)

    xrange = ohrange
    yrange = r23range

    xtitle = '12 + log (O/H) [T_{e}]'
    ytitle = 'log (R_{23})'
    
    zindicators_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, talk=dotalk, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      /left, xstyle=3, ystyle=3, /right, /top, hiipsize=1.0, position=pos[*,0]

; overlay some representative kewley grids

    hiimodel1 = interpol(hiigrids[*,1].zstrong_r23,log12ohgrid,log12ohplot,/quad)
    hiimodel2 = interpol(hiigrids[*,5].zstrong_r23,log12ohgrid,log12ohplot,/quad)
    sbmodel1 = interpol(sbgrids[*,1].zstrong_r23,log12ohgrid,log12ohplot,/quad)
    sbmodel2 = interpol(sbgrids[*,5].zstrong_r23,log12ohgrid,log12ohplot,/quad)
    
    djs_oplot, log12ohplot, hiimodel1, line=2, thick=postthick, color=djs_icolor(talkcolor)
    djs_oplot, log12ohplot, hiimodel2, line=2, thick=postthick, color=djs_icolor(talkcolor)
    djs_oplot, log12ohplot, sbmodel1, line=0, thick=postthick, color=djs_icolor(talkcolor)
    djs_oplot, log12ohplot, sbmodel2, line=0, thick=postthick, color=djs_icolor(talkcolor)
    
    xyouts, log12ohplot[0], hiimodel1[0], string(U[1],format='(F5.2)'), /data, charsize=1.5, $
      charthick=postthick, align=1.0, color=djs_icolor(talkcolor)
    xyouts, log12ohplot[0], hiimodel2[0], string(U[5],format='(F5.2)'), /data, charsize=1.5, $
      charthick=postthick, align=1.0, color=djs_icolor(talkcolor)

    legend, ['HII Region','Starburst'], line=[2,0], /right, /top, box=0, color=djs_icolor(talkcolor), $
      charsize=1.5, charthick=postthick, thick=postthick, textcolor=djs_icolor(talkcolor)

    im_openclose, postscript=postscript, /close    
    
; ###########################################################################
; electron temperature versus line ratios
; ###########################################################################

; ------------------------------------------------------------
; Te vs r23
; ------------------------------------------------------------

    psname = 'hii_te_vs_r23'
    im_openclose, pspath+psname, postscript=postscript, encapsulated=encapsulated, xsize=8.5, ysize=8.0

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.6, height=6.6, $
      xmargin=[1.5,0.4], ymargin=[0.3,1.1], xpage=8.5, ypage=8.0, $
      position=pos, /normal

    if keyword_set(dotalk) then polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, $
      color=djs_icolor('black')
    
    indx = where((hii.zstrong_r23 gt -900.0) and (hii.zt_t4363 gt -900.0),nindx)
    
    x = hii[indx].zt_t4363
    xerr = hii[indx].zt_t4363_err

    y = alog10(hii[indx].zstrong_r23)
    yerr = hii[indx].zstrong_r23_err/hii[indx].zstrong_r23/alog(10.0)

    xrange = Terange
    yrange = r23range

    xtitle = 'T_{e} [K]'
    ytitle = 'log (R_{23})'

    zindicators_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, talk=dotalk, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=3, ystyle=3, /right, /top, xtickname=Telabel, hiipsize=1.0, position=pos[*,0]
    
    im_openclose, postscript=postscript, /close    

; ------------------------------------------------------------
; Te vs [N II]/Ha
; ------------------------------------------------------------

    psname = 'hii_te_vs_niiha'
    im_openclose, pspath+psname, postscript=postscript, encapsulated=encapsulated, xsize=8.5, ysize=8.0

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.6, height=6.6, $
      xmargin=[1.5,0.4], ymargin=[0.3,1.1], xpage=8.5, ypage=8.0, $
      position=pos, /normal

    if keyword_set(dotalk) then polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, $
      color=djs_icolor('black')
    
    indx = where((hii.nii_6584_h_alpha gt -900.0) and (hii.zt_t4363 gt -900.0),nindx)
    
    x = hii[indx].zt_t4363
    xerr = hii[indx].zt_t4363_err

    y = hii[indx].nii_6584_h_alpha
    yerr = hii[indx].nii_6584_h_alpha_err

    xrange = Terange
    yrange = niiharange

    xtitle = 'T_{e} [K]'
    ytitle = 'log ([N II] \lambda6584 / H\alpha)'

    zindicators_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, talk=dotalk, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=3, ystyle=3, /right, /top, xtickname=Telabel, hiipsize=1.0, position=pos[*,0]

    im_openclose, postscript=postscript, /close    
    
; ------------------------------------------------------------
; Te vs [N II]/[O II]
; ------------------------------------------------------------

    psname = 'hii_te_vs_niioii'
    im_openclose, pspath+psname, postscript=postscript, encapsulated=encapsulated, xsize=8.5, ysize=8.0

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.6, height=6.6, $
      xmargin=[1.5,0.4], ymargin=[0.3,1.1], xpage=8.5, ypage=8.0, $
      position=pos, /normal

    if keyword_set(dotalk) then polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, $
      color=djs_icolor('black')
    
    indx = where((hii.nii_6584_oii gt -900.0) and (hii.zt_t4363 gt -900.0),nindx)
    
    x = hii[indx].zt_t4363
    xerr = hii[indx].zt_t4363_err

    y = hii[indx].nii_6584_oii
    yerr = hii[indx].nii_6584_oii_err

    xrange = Terange
    yrange = niioiirange

    xtitle = 'T_{e} [K]'
    ytitle = 'log ([N II] \lambda6584 / [O II] \lambda3727)'

    zindicators_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, talk=dotalk, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=3, ystyle=3, /right, /top, xtickname=Telabel, hiipsize=1.0, position=pos[*,0]
    
    im_openclose, postscript=postscript, /close    
    
; ------------------------------------------------------------
; Te vs o32
; ------------------------------------------------------------

    psname = 'hii_te_vs_o32'
    im_openclose, pspath+psname, postscript=postscript, encapsulated=encapsulated, xsize=8.5, ysize=8.0

    pagemaker, nx=1, ny=1, xspace=0, yspace=0.0, width=6.6, height=6.6, $
      xmargin=[1.5,0.4], ymargin=[0.3,1.1], xpage=8.5, ypage=8.0, $
      position=pos, /normal

    if keyword_set(dotalk) then polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, $
      color=djs_icolor('black')
    
    indx = where((hii.zstrong_o32 gt -900.0) and (hii.zt_t4363 gt -900.0),nindx)
    
    x = hii[indx].zt_t4363
    xerr = hii[indx].zt_t4363_err

    y = hii[indx].zstrong_o32
    yerr = hii[indx].zstrong_o32_err

    xrange = Terange
    yrange = o32range

    xtitle = 'T_{e} [K]'
    ytitle = 'log ([O III] \lambda\lambda4959,5007 / [O II] \lambda3727)'

    zindicators_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, talk=dotalk, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=3, ystyle=3, /right, /top, xtickname=Telabel, hiipsize=1.0, position=pos[*,0]
    
    im_openclose, postscript=postscript, /close    

; --------------------------------------------------    
; GENERATE TALK PLOTS
; --------------------------------------------------    

;   if keyword_set(talk) then begin
;      if (n_elements(dotalk) eq 0L) then zindicators, hii, $
;        suffix='.talk', encapsulated=0, /talk, /postscript, $
;        /dotalk, _extra=extra else return
;   endif

; --------------------------------------------------    
; GENERATE ENCAPSULATED POSTSCRIPT
; --------------------------------------------------    

;   if keyword_set(postscript) then if (n_elements(doencapsulated) eq 0L) then $
;     zindicators, hii, postscript=postscript, /encapsulated, /doencapsulated, $
;     talk=0, _extra=extra else return
    
; --------------------------------------------------    
; GENERATE THE WEB PAGE
; --------------------------------------------------    
    
    if keyword_set(postscript) then im_ps2html, htmlbase, html_path=html_path, cleanpng=0, npscols=3, _extra=extra
    
stop
    
; ###########################################################################
; the plots below should be incorporated into ABUNDANCES.PRO
; ###########################################################################

; ------------------------------------------------------------
; 12+log(O/H) ([O III]/Hb)/([N II]/Ha) vs 12+log(O/H) [Strong]
; ------------------------------------------------------------

    psname = 'hii_12oh_oiiinii_niiha_vs_12oh_strong'
    im_openclose, pspath+psname, postscript=postscript, encapsulated=encapsulated, xsize=8.5, ysize=7.8

    if keyword_set(dotalk) then polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, $
      color=djs_icolor('black')
    
    pagemaker, nx=2, ny=2, xspace=0.0, yspace=0.0, width=2.95*[1,1], height=3.35*[1,1], $
      xmargin=[1.3,1.3], ymargin=[0.2,1.3], xpage=8.5, ypage=7.8, position=pos, /normal

    xtitle = '12 + log (O/H) ([O III]/H\beta)/([N II]/H\alpha)'
    ytitle = 'Residuals'

    xrange = ohrange2
    yrange = [-0.9,0.9]

; --------------------------------------------------    
; ZKH94
; --------------------------------------------------    

    indx = where((hii.zstrong_12oh_oiiinii_pettini gt 8.1) and (hii.zstrong_12oh_zkh94 gt -900),nindx)

    x = hii[indx].zstrong_12oh_oiiinii_pettini
    xerr = hii[indx].zstrong_12oh_oiiinii_pettini_err

    y = hii[indx].zstrong_12oh_zkh94
    yerr = hii[indx].zstrong_12oh_zkh94_err

    y = x-y ; NOTE!
    
    xbig = x
    ybig = y

;   ytitle = 'log (O/H) ([O III]/H\beta)/([N II]/H\alpha) - log (O/H) [ZKH94]'
;   ytitle = '12 + log (O/H) [[N II]/[O II]]'

    stats = im_stats(ybig);,/verbose)
;   stats = im_stats(xbig-ybig,/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)
    
    rcor = r_correlate(xbig,ybig,zd=zd,probd=probd)
;   splog, 'Spearman rank test: ', rcor, zd

    zindicators_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, talk=dotalk, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=3, ystyle=3, /right, /top, position=pos[*,0], charsize=charsize_6, $
      xtickname=replicate(' ',10), noerase=keyword_set(dotalk);, xatlas=xatlas, yatlas=yatlas, $
;     xerratlas=xerratlas, yerratlas=yerratlas, xnfgs=xnfgs, ynfgs=ynfgs, $
;     xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    djs_oplot, !x.crange, [0,0], line=2, thick=postthick, color=djs_icolor(talkcolor)
;   djs_oplot, !x.crange, !y.crange, thick=postthick

    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_4, $
      charthick=postthick, textcolor=djs_icolor(talkcolor)
    legend, '(a)', /left, /top, box=0, charsize=charsize_6, charthick=postthick, $
      textcolor=djs_icolor(talkcolor)
    legend, textoidl('ZKH94'), /right, /top, box=0, charsize=charsize_6, $
      charthick=postthick, textcolor=djs_icolor(talkcolor)

; --------------------------------------------------    
; M91
; --------------------------------------------------    

    indx = where((hii.zstrong_12oh_oiiinii_pettini gt 8.1) and (hii.zstrong_12oh_m91_upper gt -900),nindx)

    x = hii[indx].zstrong_12oh_oiiinii_pettini
    xerr = hii[indx].zstrong_12oh_oiiinii_pettini_err
    
    y = hii[indx].zstrong_12oh_m91_upper
    yerr = hii[indx].zstrong_12oh_m91_upper_err

    y = x-y ; NOTE!
    
    xbig = x
    ybig = y

;   ytitle = '12 + log (O/H) M91'

    stats = im_stats(ybig);,/verbose,/no_head)
;   stats = im_stats(xbig-ybig,/verbose,/no_head)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)
    
    rcor = r_correlate(xbig,ybig,zd=zd,probd=probd)
;   splog, 'Spearman rank test: ', rcor, zd

    zindicators_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, talk=dotalk, $
      xtitle='', ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=3, ystyle=11, /right, /top, position=pos[*,1], /noerase, charsize=charsize_6, $
      xtickname=replicate(' ',10), ytickname=replicate(' ',10)
    axis, /yaxis, ystyle=3, ythick=postthick, charsize=charsize_6, charthick=postthick, $
      yrange=yrange, ytitle=ytitle, color=djs_icolor(talkcolor)
    djs_oplot, !x.crange, [0,0], line=2, thick=postthick, color=djs_icolor(talkcolor)
;   djs_oplot, !x.crange, !y.crange, thick=postthick

    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_4, $
      charthick=postthick, textcolor=djs_icolor(talkcolor)
    legend, '(b)', /left, /top, box=0, charsize=charsize_6, charthick=postthick, $
      textcolor=djs_icolor(talkcolor)
    legend, textoidl('M91'), /right, /top, box=0, charsize=charsize_6, $
      charthick=postthick, textcolor=djs_icolor(talkcolor)
    
; --------------------------------------------------    
; KD02-Combined
; --------------------------------------------------    

    indx = where((hii.zstrong_12oh_oiiinii_pettini gt 8.1) and (hii.zstrong_12oh_kd02_combined gt -900.0),nindx)
 
    x = hii[indx].zstrong_12oh_oiiinii_pettini
    xerr = hii[indx].zstrong_12oh_oiiinii_pettini_err

    y = hii[indx].zstrong_12oh_kd02_combined
    yerr = hii[indx].zstrong_12oh_kd02_combined_err

    y = x-y ; NOTE!

    xbig = x
    ybig = y

;   ytitle = '12 + log (O/H) [KD02]'

    stats = im_stats(ybig);,/verbose,/no_head)
;   stats = im_stats(xbig-ybig,/verbose,/no_head)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)
    
    rcor = r_correlate(xbig,ybig,zd=zd,probd=probd)
;   splog, 'Spearman rank test: ', rcor, zd

    zindicators_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, talk=dotalk, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, charsize=charsize_6, $
      xstyle=3, ystyle=3, /right, /top, position=pos[*,2], /noerase;, xnfgs=xnfgs, ynfgs=ynfgs, $
;     xerrnfgs=xerrnfgs, yerrnfgs=yerrnfgs
    djs_oplot, !x.crange, [0,0], line=2, thick=postthick, color=djs_icolor(talkcolor)
;   djs_oplot, !x.crange, !y.crange, thick=postthick
 
    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_4, $
      charthick=postthick, textcolor=djs_icolor(talkcolor)
    legend, '(c)', /left, /top, box=0, charsize=charsize_6, charthick=postthick, $
      textcolor=djs_icolor(talkcolor)
    legend, textoidl('KD02'), /right, /top, box=0, charsize=charsize_6, $
      charthick=postthick, textcolor=djs_icolor(talkcolor)

; --------------------------------------------------    
; T04
; --------------------------------------------------    

    indx = where((hii.zstrong_12oh_oiiinii_pettini gt 8.1) and (hii.zstrong_12oh_t04 gt -900),nindx)

    x = hii[indx].zstrong_12oh_oiiinii_pettini
    xerr = hii[indx].zstrong_12oh_oiiinii_pettini_err
    
    y = hii[indx].zstrong_12oh_t04
    yerr = hii[indx].zstrong_12oh_t04_err

    y = x-y ; NOTE!
    
    xbig = x
    ybig = y

;   ytitle = '12 + log (O/H) [T04]'

    stats = im_stats(ybig);,/verbose,/no_head)
;   stats = im_stats(xbig-ybig,/verbose,/no_head)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)
    
    rcor = r_correlate(xbig,ybig,zd=zd,probd=probd)
;   splog, 'Spearman rank test: ', rcor, zd

    zindicators_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, talk=dotalk, $
      xtitle='', ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=3, ystyle=11, /right, /top, position=pos[*,3], /noerase, charsize=charsize_6, $
      ytickname=replicate(' ',10)
    axis, /yaxis, ystyle=3, ythick=postthick, charsize=charsize_6, charthick=postthick, $
      yrange=yrange, ytitle=ytitle, color=djs_icolor(talkcolor)
    djs_oplot, !x.crange, [0,0], line=2, thick=postthick, color=djs_icolor(talkcolor)

    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_4, $
      charthick=postthick, textcolor=djs_icolor(talkcolor)
    legend, '(d)', /left, /top, box=0, charsize=charsize_6, charthick=postthick, $
      textcolor=djs_icolor(talkcolor)
    legend, textoidl('T04'), /right, /top, box=0, charsize=charsize_6, $
      charthick=postthick, textcolor=djs_icolor(talkcolor)

; title    
    
    xyouts, 0.5*(pos[2,1]-pos[0,0])+pos[0,0], pos[1,2]*0.2, textoidl(xtitle), align=0.5, $
      charsize=charsize_6, charthick=postthick, /normal, color=djs_icolor(talkcolor)

    im_openclose, postscript=postscript, /close    
    
; ------------------------------------------------------------
; inter-compare empirical and strong-line calibrations
; ------------------------------------------------------------

    psname = 'hii_12oh_strong_vs_12oh_empirical'
    im_openclose, pspath+psname, postscript=postscript, encapsulated=encapsulated

    pagemaker, nx=3, ny=3, position=pos, /normal, xspace=0.0, $
      xmargin=[1.3,0.2], ymargin=[0.2,1.1], yspace=0.0

    if keyword_set(dotalk) then polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, $
      color=djs_icolor('black')
    
    xrange = ohrange3
    yrange = xrange
    charsize_6 = 1.2

; --------------------------------------------------    
; Panel 1 - oiiinii+niiha versus KD02-Combined
; --------------------------------------------------    

    indx = where((hii.zstrong_12oh_oiiinii_pettini gt -900.0) and (hii.zstrong_12oh_kd02_combined gt -900.0),nindx)
 
    x = hii[indx].zstrong_12oh_oiiinii_pettini
    xerr = hii[indx].zstrong_12oh_oiiinii_pettini_err

    y = hii[indx].zstrong_12oh_kd02_combined
    yerr = hii[indx].zstrong_12oh_kd02_combined_err

    xtitle = '12 + log (O/H) ([O III]/H\beta)/([N II]/H\alpha)'
    ytitle = '12 + log (O/H) KD02'

    stats = im_stats(x-y,/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)
    
    zindicators_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, talk=dotalk, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, charsize=charsize_6, $
      xstyle=3, ystyle=3, /right, /top, position=pos[*,0], xtickname=replicate(' ',10), $
      noerase=keyword_set(dotalk)
    djs_oplot, !x.crange, !y.crange, thick=postthick, color=djs_icolor(talkcolor)
 
    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_6, charthick=postthick, $
      textcolor=djs_icolor(talkcolor)
    legend, '(a)', /left, /top, box=0, charsize=charsize_6, charthick=postthick, textcolor=djs_icolor(talkcolor)
    
; --------------------------------------------------    
; Panel 2 - niihaO2 versus KD02-Combined
; --------------------------------------------------    

    indx = where((hii.zstrong_12oh_kd02_niioii gt -900.0) and (hii.zstrong_12oh_kd02_combined gt -900.0),nindx)

    x = hii[indx].zstrong_12oh_kd02_niioii
    xerr = hii[indx].zstrong_12oh_kd02_niioii_err

    y = hii[indx].zstrong_12oh_kd02_combined
    yerr = hii[indx].zstrong_12oh_kd02_combined_err

    xtitle = '12 + log (O/H) [N II]/[O II]'
    ytitle = '12 + log (O/H) KD02'

    stats = im_stats(x-y,/verbose,/no_head)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)
    
    zindicators_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, talk=dotalk, $
      xtitle='', ytitle='', xrange=xrange, yrange=yrange, legendtype=0, charsize=charsize_6, $
      xstyle=3, ystyle=3, /right, /top, position=pos[*,1], /noerase, ytickname=replicate(' ',10), $
      xtickname=replicate(' ',10)
    djs_oplot, !x.crange, !y.crange, thick=postthick, color=djs_icolor(talkcolor)

    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_6, charthick=postthick, $
      textcolor=djs_icolor(talkcolor)
    legend, '(b)', /left, /top, box=0, charsize=charsize_6, charthick=postthick, textcolor=djs_icolor(talkcolor)

; --------------------------------------------------    
; Panel 3 - M91 versus KD02-Combined
; --------------------------------------------------    
 
    indx = where((hii.zstrong_12oh_m91_upper gt -900.0) and (hii.zstrong_12oh_m91_lower gt -900.0) and $
      (hii.zstrong_12oh_kd02_combined gt -900.0),nindx)

    y = hii[indx].zstrong_12oh_kd02_combined
    yerr = hii[indx].zstrong_12oh_kd02_combined_err

    x = y*0.0
    xerr = yerr*0.0
    
    for idiv = 0L, ndiv-1L do begin

       up = where((hii[indx].zstrong_12oh_kd02_combined gt div[idiv]),nup,comp=lo,ncomp=nlo)
       if (nup ne 0L) then begin
          x[up] = hii[indx[up]].zstrong_12oh_m91_upper
          xerr[up] = hii[indx[up]].zstrong_12oh_m91_upper_err
       endif

       if (nlo ne 0L) then begin
          x[lo] = hii[indx[lo]].zstrong_12oh_m91_lower
          xerr[lo] = hii[indx[lo]].zstrong_12oh_m91_lower_err
       endif

       resid[idiv] = stddev(x-y)
;      resid[idiv] = djsig(x-y,sigrej=3.0)
;      plot, x, x-y, ps=4, xrange=xrange, yrange=[-1,1]
;      print, div[idiv], resid[idiv]
;      cc = get_kbrd(1)
       
    endfor

;   mindiv = find_nminima(resid,div,nfind=1,width=0.2)
    minresid = min(resid,minindx)
    mindiv = div[minindx]
    splog, 'M91 vs KD02-Combined: '+string(mindiv,format='(F4.2)')

    up = where((hii[indx].zstrong_12oh_kd02_combined gt mindiv),nup)
    if (nup ne 0L) then begin
       x = hii[indx[up]].zstrong_12oh_m91_upper
       xerr = hii[indx[up]].zstrong_12oh_m91_upper_err

       y = hii[indx[up]].zstrong_12oh_kd02_combined
       yerr = hii[indx[up]].zstrong_12oh_kd02_combined_err
    endif

    lo = where((hii[indx].zstrong_12oh_kd02_combined lt mindiv),nlo)
    if (nlo ne 0L) then begin
       x = [x,hii[indx[lo]].zstrong_12oh_m91_lower]
       xerr = [xerr,hii[indx[lo]].zstrong_12oh_m91_lower_err]

       y = [y,hii[indx[lo]].zstrong_12oh_kd02_combined]
       yerr = [yerr,hii[indx[lo]].zstrong_12oh_kd02_combined_err]
    endif

    xtitle = '12 + log (O/H) M91'
    ytitle = '12 + log (O/H) KD02'

    stats = im_stats(x-y,/verbose,/no_head)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)
    
    zindicators_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, talk=dotalk, $
      xtitle=xtitle, ytitle='', xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=3, ystyle=3, /right, /top, position=pos[*,2], /noerase, charsize=charsize_6, $
      ytickname=replicate(' ',10)
    djs_oplot, !x.crange, !y.crange, thick=postthick, color=djs_icolor(talkcolor)
 
    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_6, charthick=postthick, $
      textcolor=djs_icolor(talkcolor)
    legend, '(c)', /left, /top, box=0, charsize=charsize_6, charthick=postthick, textcolor=djs_icolor(talkcolor)

; --------------------------------------------------    
; Panel 4 - oiiinii+niiha versus M91
; --------------------------------------------------    

    indx = where((hii.zstrong_12oh_oiiinii_pettini gt -900) and (hii.zstrong_12oh_m91_upper gt -900.0) and $
      (hii.zstrong_12oh_m91_lower gt -900.0),nindx)

    x = hii[indx].zstrong_12oh_oiiinii_pettini
    xerr = hii[indx].zstrong_12oh_oiiinii_pettini_err
    
    y = x*0.0
    yerr = xerr*0.0
    
    for idiv = 0L, ndiv-1L do begin

       up = where((x gt div[idiv]),nup,comp=lo,ncomp=nlo)
       if (nup ne 0L) then begin
          y[up] = hii[indx[up]].zstrong_12oh_m91_upper
          yerr[up] = hii[indx[up]].zstrong_12oh_m91_upper_err
       endif

       if (nlo ne 0L) then begin
          y[lo] = hii[indx[lo]].zstrong_12oh_m91_lower
          yerr[lo] = hii[indx[lo]].zstrong_12oh_m91_lower_err
       endif

       resid[idiv] = stddev(x-y)
;      resid[idiv] = djsig(x-y,sigrej=3.0)
;      plot, x, x-y, ps=4, xrange=xrange, yrange=[-1,1]
;      print, div[idiv], resid[idiv]
;      cc = get_kbrd(1)
       
    endfor

;   mindiv = find_nminima(resid,div,nfind=1,width=0.2)
    minresid = min(resid,minindx)
    mindiv = div[minindx]
    splog, 'oiiinii+niiha vs M91: '+string(mindiv,format='(F4.2)')

    up = where((hii[indx].zstrong_12oh_oiiinii_pettini gt mindiv),nup)
    if (nup ne 0L) then begin
       x = hii[indx[up]].zstrong_12oh_oiiinii_pettini
       xerr = hii[indx[up]].zstrong_12oh_oiiinii_pettini_err

       y = hii[indx[up]].zstrong_12oh_m91_upper
       yerr = hii[indx[up]].zstrong_12oh_m91_upper_err
    endif

    lo = where((hii[indx].zstrong_12oh_oiiinii_pettini lt mindiv),nlo)
    if (nlo ne 0L) then begin
       x = [x,hii[indx[lo]].zstrong_12oh_oiiinii_pettini]
       xerr = [xerr,hii[indx[lo]].zstrong_12oh_oiiinii_pettini_err]

       y = [y,hii[indx[lo]].zstrong_12oh_m91_lower]
       yerr = [yerr,hii[indx[lo]].zstrong_12oh_m91_lower_err]
    endif

    xtitle = '12 + log (O/H) ([O III]/H\beta)/([N II]/H\alpha)'
    ytitle = '12 + log (O/H) M91'

    stats = im_stats(x-y,/verbose,/no_head)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)
    
    zindicators_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, talk=dotalk, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=3, ystyle=3, /right, /top, position=pos[*,3], /noerase, charsize=charsize_6, $
      xtickname=replicate(' ',10)
    djs_oplot, !x.crange, !y.crange, thick=postthick, color=djs_icolor(talkcolor)

    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_6, charthick=postthick, $
      textcolor=djs_icolor(talkcolor)
    legend, '(d)', /left, /top, box=0, charsize=charsize_6, charthick=postthick, textcolor=djs_icolor(talkcolor)
    
; --------------------------------------------------    
; Panel 5 - niihaO2 versus M91
; --------------------------------------------------    

    indx = where((hii.zstrong_12oh_kd02_niioii gt -900.0) and (hii.zstrong_12oh_m91_upper gt -900.0) and $
      (hii.zstrong_12oh_m91_lower gt -900.0),nindx)
    
    x = hii[indx].zstrong_12oh_kd02_niioii
    xerr = hii[indx].zstrong_12oh_kd02_niioii_err
    
    y = x*0.0
    yerr = xerr*0.0
    
    for idiv = 0L, ndiv-1L do begin

       up = where((x gt div[idiv]),nup,comp=lo,ncomp=nlo)
       if (nup ne 0L) then begin
          y[up] = hii[indx[up]].zstrong_12oh_m91_upper
          yerr[up] = hii[indx[up]].zstrong_12oh_m91_upper_err
       endif

       if (nlo ne 0L) then begin
          y[lo] = hii[indx[lo]].zstrong_12oh_m91_lower
          yerr[lo] = hii[indx[lo]].zstrong_12oh_m91_lower_err
       endif

       resid[idiv] = stddev(x-y)
;      resid[idiv] = djsig(x-y,sigrej=3.0)
;      plot, x, x-y, ps=4, xrange=xrange, yrange=[-1,1]
;      print, div[idiv], resid[idiv]
;      cc = get_kbrd(1)
       
    endfor

;   mindiv = find_nminima(resid,div,nfind=1,width=0.2)
    minresid = min(resid,minindx)
    mindiv = div[minindx]
    splog, 'niihaO2 vs M91: '+string(mindiv,format='(F4.2)')

    up = where((hii[indx].zstrong_12oh_kd02_niioii gt mindiv),nup)
    if (nup ne 0L) then begin
       x = hii[indx[up]].zstrong_12oh_kd02_niioii
       xerr = hii[indx[up]].zstrong_12oh_kd02_niioii_err

       y = hii[indx[up]].zstrong_12oh_m91_upper
       yerr = hii[indx[up]].zstrong_12oh_m91_upper_err
    endif

    lo = where((hii[indx].zstrong_12oh_kd02_niioii lt mindiv),nlo)
    if (nlo ne 0L) then begin
       x = [x,hii[indx[lo]].zstrong_12oh_kd02_niioii]
       xerr = [xerr,hii[indx[lo]].zstrong_12oh_kd02_niioii_err]

       y = [y,hii[indx[lo]].zstrong_12oh_m91_lower]
       yerr = [yerr,hii[indx[lo]].zstrong_12oh_m91_lower_err]
    endif

    xtitle = '12 + log (O/H) [N II]/[O II]'
    ytitle = '12 + log (O/H) M91'

    stats = im_stats(x-y,/verbose,/no_head)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)
    
    zindicators_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, talk=dotalk, $
      xtitle=xtitle, ytitle='', xrange=xrange, yrange=yrange, legendtype=0, charsize=charsize_6, $
      xstyle=3, ystyle=3, /right, /top, position=pos[*,4], /noerase, $
      ytickname=replicate(' ',10)
    djs_oplot, !x.crange, !y.crange, thick=postthick, color=djs_icolor(talkcolor)

    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_6, charthick=postthick, $
      textcolor=djs_icolor(talkcolor)
    legend, '(e)', /left, /top, box=0, charsize=charsize_6, charthick=postthick, textcolor=djs_icolor(talkcolor)
    
; --------------------------------------------------    
; Panel 6 - No Data
; --------------------------------------------------    

; --------------------------------------------------    
; Panel 7 - oiiinii+niiha versus niihaO2
; --------------------------------------------------    

    indx = where((hii.zstrong_12oh_oiiinii_pettini gt -900) and (hii.zstrong_12oh_kd02_niioii gt -900.0),nindx)

    x = hii[indx].zstrong_12oh_oiiinii_pettini
    xerr = hii[indx].zstrong_12oh_oiiinii_pettini_err

    y = hii[indx].zstrong_12oh_kd02_niioii
    yerr = hii[indx].zstrong_12oh_kd02_niioii_err

    xtitle = '12 + log (O/H) ([O III]/H\beta)/([N II]/H\alpha)'
    ytitle = '12 + log (O/H) [N II]/[O II]'

    stats = im_stats(x-y,/verbose,/no_head)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)
    
    zindicators_lineplot, x, y, xerr, yerr, plottype=1, postscript=postscript, talk=dotalk, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=3, ystyle=3, /right, /top, position=pos[*,6], /noerase, charsize=charsize_6
    djs_oplot, !x.crange, !y.crange, thick=postthick, color=djs_icolor(talkcolor)

    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_6, charthick=postthick, $
      textcolor=djs_icolor(talkcolor)
    legend, '(f)', /left, /top, box=0, charsize=charsize_6, charthick=postthick, textcolor=djs_icolor(talkcolor)
    
; --------------------------------------------------    
; Panel 8 - No Data
; --------------------------------------------------    

; --------------------------------------------------    
; Panel 9 - No Data
; --------------------------------------------------    

    im_openclose, postscript=postscript, /close    

; ###########################################################################
; the plots above should be incorporated into ABUNDANCES.PRO
; ###########################################################################

; ###########################################################################
; the next two plots need more symbols defined in IM_SYMBOLS()!
; ###########################################################################

; ###########################################################################
; the following plots need to be fixed or are outdated
; ###########################################################################

; ------------------------------------------------------------
; 12+log(O/H) Te vs 12+log(O/H) Empirical - Residuals
; ------------------------------------------------------------

    psname = 'hii_12oh_te_vs_12oh_empirical_residuals'
    im_openclose, pspath+psname, postscript=postscript, encapsulated=encapsulated, xsize=8.5, ysize=10.5

    pagemaker, nx=1, ny=3, height=3.0*[1,1,1], width=6.8, xmargin=[1.3,0.4], $
      ymargin=[0.4,1.1], xspace=0, yspace=0, xpage=8.5, ypage=10.5, $
      position=pos, /normal

    if keyword_set(dotalk) then polyfill, [1,1,0,0,1], [1,0,0,1,1], /normal, $
      color=djs_icolor('black')

    ytitle = '\Delta[log(O/H)]'
    xtitle = ohtitle+' T_{e}'

    xrange = ohrange
    yrange = [-1.2,1.2]
    
; --------------------------------------------------    
; ([O III]/Hb)/([N II]/Ha)
; --------------------------------------------------    

    indx = where((hii.zstrong_12oh_oiiinii_moustakas gt -900) and (hii.zt_log12oh_te gt -900.0),nindx)

    x = hii[indx].zt_log12oh_te
    xerr = hii[indx].zt_log12oh_te_err

    y1 = hii[indx].zstrong_12oh_oiiinii_moustakas
    y1err = hii[indx].zstrong_12oh_oiiinii_moustakas_err

    y = x-y1
    yerr = sqrt(x^2 + y1err^2)
    
    stats = im_stats(y,/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)
    
    zindicators_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, talk=dotalk, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=3, ystyle=3, /right, /top, position=pos[*,0], charsize=charsize_6, $
      xtickname=replicate(' ',10), hiipsize=0.8
    djs_oplot, !x.crange, [0,0], thick=postthick, color=djs_icolor(talkcolor), line=2

    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_5, $
      charthick=postthick, textcolor=djs_icolor(talkcolor)
    legend, textoidl('(a) ([O III]/H\beta)/([N II]/H\alpha)'), /left, /top, box=0, $
      charsize=charsize_5, charthick=postthick, textcolor=djs_icolor(talkcolor)

; --------------------------------------------------    
; [N II]/Ha
; --------------------------------------------------    

    indx = where((hii.zstrong_12oh_niiha_moustakas gt -900) and (hii.zt_log12oh_te gt -900.0),nindx)

    x = hii[indx].zt_log12oh_te
    xerr = hii[indx].zt_log12oh_te_err

    y1 = hii[indx].zstrong_12oh_niiha_moustakas
    y1err = hii[indx].zstrong_12oh_niiha_moustakas_err

    y = x-y1
    yerr = sqrt(x^2 + y1err^2)
    
    stats = im_stats(y,/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)
    
    zindicators_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, talk=dotalk, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=3, ystyle=3, /right, /top, position=pos[*,1], charsize=charsize_6, $
      xtickname=replicate(' ',10), hiipsize=0.8, /noerase
    djs_oplot, !x.crange, [0,0], thick=postthick, color=djs_icolor(talkcolor), line=2

    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_5, $
      charthick=postthick, textcolor=djs_icolor(talkcolor)
    legend, textoidl('(b) [N II]/H\alpha'), /left, /top, box=0, $
      charsize=charsize_5, charthick=postthick, textcolor=djs_icolor(talkcolor)

; --------------------------------------------------    
; P-method
; --------------------------------------------------    

    indx = where((hii.zt_log12oh_te gt -900.0) and (hii.zstrong_12oh_p01_upper gt -900.0) and $
      (hii.zstrong_12oh_p01_lower gt -900.0),nindx)

    up = where((hii[indx].zt_log12oh_te gt 8.2),nup)
    if (nup ne 0L) then begin
       x = hii[indx[up]].zt_log12oh_te
       xerr = hii[indx[up]].zt_log12oh_te_err

       y1 = hii[indx[up]].zstrong_12oh_p01_upper
       y1err = hii[indx[up]].zstrong_12oh_p01_upper_err
    endif

    lo = where((hii[indx].zt_log12oh_te lt 8.0),nlo)
    if (nlo ne 0L) then begin
       x = [x,hii[indx[lo]].zt_log12oh_te]
       xerr = [xerr,hii[indx[lo]].zt_log12oh_te_err]

       y1 = [y1,hii[indx[lo]].zstrong_12oh_p01_lower]
       y1err = [y1err,hii[indx[lo]].zstrong_12oh_p01_lower_err]
    endif

    y = x-y1
    yerr = sqrt(x^2 + y1err^2)
    
    stats = im_stats(y,/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)
    
    zindicators_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, talk=dotalk, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=3, ystyle=3, /right, /top, position=pos[*,2], charsize=charsize_6, $
      hiipsize=0.8, /noerase
    djs_oplot, !x.crange, [0,0], thick=postthick, color=djs_icolor(talkcolor), line=2

    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_5, $
      charthick=postthick, textcolor=djs_icolor(talkcolor)
    legend, textoidl('(c) P-method'), /left, /top, box=0, $
      charsize=charsize_5, charthick=postthick, textcolor=djs_icolor(talkcolor)

    im_openclose, postscript=postscript, /close    
    
; ------------------------------------------------------------
; [N II]/Ha vs [O III]/Hb - instantaneous
; ------------------------------------------------------------

    psname = 'hii_models_inst_niiha_vs_oiiihb'
    im_openclose, pspath+psname, postscript=postscript

    plot_kewley_grids, plotnumber=1, model=8, label=1, /hii, postscript=postscript
    legend, 'SB99 - instantaneous burst', /left, /bottom, box=0, charsize=1.5, charthick=postthick
    
    im_openclose, postscript=postscript, /close
    
; ------------------------------------------------------------
; [N II]/Ha vs [O III]/Hb - continuous
; ------------------------------------------------------------

    psname = 'hii_models_cont_niiha_vs_oiiihb'
    im_openclose, pspath+psname, postscript=postscript

    plot_kewley_grids, plotnumber=1, model=3, label=1, /hii, postscript=postscript
    legend, 'SB99 - continous burst', /left, /bottom, box=0, charsize=1.5, charthick=postthick
    
    im_openclose, postscript=postscript, /close
    
; ------------------------------------------------------------
; r23 vs [N II]/[O II] - instantaneous
; ------------------------------------------------------------

    psname = 'hii_models_inst_r23_vs_niioii'
    im_openclose, pspath+psname, postscript=postscript

    plot_kewley_grids, plotnumber=8, model=8, label=1, /hii, postscript=postscript
    legend, 'SB99 - instantaneous burst', /left, /bottom, box=0, charsize=1.5, charthick=postthick
    
    im_openclose, postscript=postscript, /close
    
; ------------------------------------------------------------
; r23 vs [N II]/[O II] - continuous
; ------------------------------------------------------------

    psname = 'hii_models_cont_r23_vs_niioii'
    im_openclose, pspath+psname, postscript=postscript

    plot_kewley_grids, plotnumber=8, model=3, label=1, /hii, postscript=postscript
    legend, 'SB99 - continous burst', /left, /bottom, box=0, charsize=1.5, charthick=postthick
    
    im_openclose, postscript=postscript, /close
    
; ------------------------------------------------------------
; r23 vs [O III]/[O II] - instantaneous
; ------------------------------------------------------------

    psname = 'hii_models_inst_r23_vs_oiiioii'
    im_openclose, pspath+psname, postscript=postscript

    plot_kewley_grids, plotnumber=10, model=8, label=1, /hii, postscript=postscript
    legend, 'SB99 - instantaneous burst', /left, /bottom, box=0, charsize=1.5, charthick=postthick
    
    im_openclose, postscript=postscript, /close
    
; ------------------------------------------------------------
; r23 vs [O III]/[O II] - continuous
; ------------------------------------------------------------

    psname = 'hii_models_cont_r23_vs_oiiioii'
    im_openclose, pspath+psname, postscript=postscript

    plot_kewley_grids, plotnumber=10, model=3, label=1, /hii, postscript=postscript
    legend, 'SB99 - continuous burst', /left, /bottom, box=0, charsize=1.5, charthick=postthick
    
    im_openclose, postscript=postscript, /close
    
; ------------------------------------------------------------
; [N II]/Ha vs [O III]/[O II] - instantaneous
; ------------------------------------------------------------

    psname = 'hii_models_inst_niiha_vs_oiiioii'
    im_openclose, pspath+psname, postscript=postscript

    plot_kewley_grids, plotnumber=23, model=8, label=1, /hii, postscript=postscript
    legend, 'SB99 - instantaneous burst', /left, /bottom, box=0, charsize=1.5, charthick=postthick
    
    im_openclose, postscript=postscript, /close
    
; ------------------------------------------------------------
; [N II]/Ha vs [O III]/[O II] - continuous
; ------------------------------------------------------------

    psname = 'hii_models_cont_niiha_vs_oiiioii'
    im_openclose, pspath+psname, postscript=postscript

    plot_kewley_grids, plotnumber=23, model=3, label=1, /hii, postscript=postscript
    legend, 'SB99 - continuous burst', /left, /bottom, box=0, charsize=1.5, charthick=postthick
    
    im_openclose, postscript=postscript, /close

; ------------------------------------------------------------
; 12 + log (O/H) [Te] vs 12 + log (O/H) [P01] - Plot II
; ------------------------------------------------------------

    psname = '12oh_P01_vs_12oh_residuals'
    im_openclose, pspath+psname, postscript=postscript

    lower = where((hii.zstrong_12oh_P01_lower gt -900.0) and (hii.zt_log12oh_te gt -900.0) and $
      (hii.zt_log12oh_te le 8.0),nlower)
    upper = where((hii.zstrong_12oh_P01_upper gt -900.0) and (hii.zt_log12oh_te gt -900.0) and $
      (hii.zt_log12oh_te ge 8.2),nupper)

    indx = cmset_op(lower,'OR',upper)

; -----
    Plo = where(hii[lower].zstrong_P lt 0.4,comp=Phi)

    xlower_Plo = hii[lower[Plo]].zstrong_12oh_P01_lower
    xlowererr_Plo = hii[lower[Plo]].zstrong_12oh_P01_lower_err

    xlower_Phi = hii[lower[Phi]].zstrong_12oh_P01_lower
    xlowererr_Phi = hii[lower[Phi]].zstrong_12oh_P01_lower_err

    ylower_Plo = hii[lower[Plo]].zstrong_12oh_P01_lower - hii[lower[Plo]].zt_log12oh_te
    ylowererr_Plo = sqrt(hii[lower[Plo]].zstrong_12oh_P01_lower_err^2 + hii[lower[Plo]].zt_log12oh_te_err^2)

    ylower_Phi = hii[lower[Phi]].zstrong_12oh_P01_lower - hii[lower[Phi]].zt_log12oh_te
    ylowererr_Phi = sqrt(hii[lower[Phi]].zstrong_12oh_P01_lower_err^2 + hii[lower[Phi]].zt_log12oh_te_err^2)
; -----

    Plo = where(hii[upper].zstrong_P lt 0.4,comp=Phi)

    xupper_Plo = hii[upper[Plo]].zstrong_12oh_P01_upper
    xuppererr_Plo = hii[upper[Plo]].zstrong_12oh_P01_upper_err

    xupper_Phi = hii[upper[Phi]].zstrong_12oh_P01_upper
    xuppererr_Phi = hii[upper[Phi]].zstrong_12oh_P01_upper_err

    yupper_Plo = hii[upper[Plo]].zstrong_12oh_P01_upper - hii[upper[Plo]].zt_log12oh_te
    yuppererr_Plo = sqrt(hii[upper[Plo]].zstrong_12oh_P01_upper_err^2 + hii[upper[Plo]].zt_log12oh_te_err^2)

    yupper_Phi = hii[upper[Phi]].zstrong_12oh_P01_upper - hii[upper[Phi]].zt_log12oh_te
    yuppererr_Phi = sqrt(hii[upper[Phi]].zstrong_12oh_P01_upper_err^2 + hii[upper[Phi]].zt_log12oh_te_err^2)
; -----

    xrange = ohrange
;   yrange = [-0.75,0.75]
    yrange = [-1.2,1.2]

    xtitle = '12 + log (O/H) [P01}]'
    ytitle = 'log (O/H) [P01] - log (O/H) T_{e}'

    djs_plot, [0], [0], /nodata, xtitle=xtitle, ytitle=ytitle, $
      xrange=xrange, yrange=yrange, xstyle=3, ystyle=3, $
      xthick=postthick, ythick=postthick, charsize=2.0, charthick=postthick
    djs_oplot, !x.crange, [0,0], thick=postthick

    plotsym, 0, 1, /fill
    djs_oplot, [xlower_Phi,xupper_Phi], [ylower_Phi,yupper_Phi], ps=8, color='red'

    plotsym, 8, 1, /fill
    djs_oplot, [xlower_Plo,xupper_Plo], [ylower_Plo,yupper_Plo], ps=8, color='green'
    
    im_legend, ['P > 0.4','P < 0.4'], /right, /top, box=0, charsize=2.0, $
      charthick=postthick, psym=[108,106], fill=[1,1], symsize=1.4, $
      color=djs_icolor(['red','green'])

; print some statistics
    
    ylower_big = [ylower_Plo,ylower_Phi]
    yupper_big = [yupper_Plo,yupper_Phi]

    splog, 'LOWER, UPPER, COMBINED statistics:'
    stats = im_stats(ylower_big,/verbose)
    stats = im_stats(yupper_big,/verbose,/no_head)
    stats = im_stats([ylower_big,yupper_big],/verbose,/no_head)

    im_openclose, postscript=postscript, /close    


; ------------------------------------------------------------
; 2D calibration of abundance using niiha and [O III]/[O II]
; ------------------------------------------------------------
;
;   psname = 'test'
;   im_openclose, pspath+psname, postscript=postscript
;
;   indx = where((hii.zstrong_niiha gt -900.0) and (hii.oiii_5007_oii gt -900.0) and $
;     (hii.zt_log12oh_te gt -900.0),nindx)
;    
;   x = hii[indx].zstrong_niiha
;   xerr = hii[indx].zstrong_niiha_err
;   x2d = hii[indx].zstrong_niiha # (replicate(1.0,nindx))
;
;   y = hii[indx].oiii_5007_oii
;   yerr = hii[indx].oiii_5007_oii_err
;   y2d = hii[indx].oiii_5007_oii # (replicate(1.0,nindx))
;
;   ndegree = 4L
;   result = sfit([ [x], [y] ],ndegree,kx=kx)
;   surface, result, x2d, y2d
;   
;   z = hii[indx].zt_log12oh_te
;   zerr = hii[indx].zt_log12oh_te_err
;
;   z2d = z # (dblarr(nindx)+1)
;   
;   xrange = [-0.5,1.2]
;   yrange = [-1.5,1.5]
;   zrange = ohrange
;
;   xtitle = textoidl('log (R_{23})')
;   ytitle = textoidl('log ([O III] / [O II])')
;   ztitle = textoidl('12 + log (O/H)')
;
;   plotsym, 0, 1.0, /fill
;   plot_3dbox, x, y, z, ps=8, xsty=3, ysty=3, zsty=3, charsize=3.0, $
;     xrange=xrange, yrange=yrange, zrange=zrange, xtitle=xtitle, $
;     ytitle=ytitle, ztitle=ztitle, charthick=postthick, xthick=postthick, $
;     ythick=postthick, zthick=postthick, gridstyle=1
;
;   im_openclose, postscript=postscript, /close    
    
; ------------------------------------------------------------
; 
; ------------------------------------------------------------
 
;   psname = ''
;   im_openclose, pspath+psname, postscript=postscript
;
;   indx = where((hii.zstrong_r23 gt -900.0) and (hii.oiii_5007_oii gt -900.0) and $
;     ((hii.zstrong_12oh_n2 gt -900.0) or (hii.zt_log12oh_te gt -900.0)),nindx)
;    
;   x = alog10(hii[indx].zstrong_r23)
;   xerr = hii[indx].zstrong_r23_err/alog10(hii[indx].zstrong_r23)/alog(10.0)
;
;   y = hii[indx].oiii_5007_oii
;   yerr = hii[indx].oiii_5007_oii_err
;
;   z1 = transpose([ [hii[indx].zstrong_12oh_n2], [hii[indx].zt_log12oh_te] ])
;   z1err = transpose([ [hii[indx].zstrong_12oh_n2_err], [hii[indx].zt_log12oh_te_err] ])
;
;   z = x*0.0
;   zerr = xerr*0.0
;
;   for i = 0L, nindx-1L do if (z1[1,i] gt 0.0) then z[i] = z1[1,i] else z[i] = z1[0,i]
;   for i = 0L, nindx-1L do if (z1[1,i] gt 0.0) then zerr[i] = z1err[1,i] else zerr[i] = z1err[0,i]
;
;   z2d = z # (dblarr(nindx)+1)
;   
;   xrange = [-0.5,1.2]
;   yrange = [-1.5,1.5]
;   zrange = ohrange
;
;   xtitle = textoidl('log (R_{23})')
;   ytitle = textoidl('log ([O III] / [O II])')
;   ztitle = textoidl('12 + log (O/H)')
;
;   plotsym, 0, 1.0, /fill
;   plot_3dbox, x, y, z, ps=8, xsty=3, ysty=3, zsty=3, charsize=3.0, $
;     xrange=xrange, yrange=yrange, zrange=zrange, xtitle=xtitle, $
;     ytitle=ytitle, ztitle=ztitle, charthick=postthick, xthick=postthick, $
;     ythick=postthick, zthick=postthick, gridstyle=1
;
;   im_openclose, postscript=postscript, /close    


; ###########################################################################
; the following plots are outdated
; ###########################################################################

    
    
; ------------------------------------------------------------
; Various Strong-line Abundances vs Te-based abundances - Plot I
; ------------------------------------------------------------

    psname = '12oh_strong_vs_12oh_residuals'
    im_openclose, pspath+psname, postscript=postscript

    pagemaker, nx=2, ny=2, position=pos, /normal, xspace=0.0, $
      xmargin=[1.6,0.1], ymargin=[0.5,1.2], yspace=0.0

    ytitle = 'log (O/H) [Strong-line] - log (O/H) T_{e}'
    xtitle = '12 + log (O/H) [Strong-line]'

;   xrange = [7.1,9.4]
    xrange = ohrange
    yrange = [-1.2,1.2]
    charsize_6 = 1.8

; --------------------------------------------------    
; Charlot & Longhetti (2001) - Case A
; --------------------------------------------------    

    indx = where(hii.zstrong_12oh_CL01A gt -900.0 and hii.zt_log12oh_te gt -900.0,nindx)

    x = hii[indx].zstrong_12oh_CL01A
;   xerr = hii[indx].zstrong_12oh_CL01A_err
    xerr = x*0.0
    
    y = hii[indx].zstrong_12oh_CL01A-hii[indx].zt_log12oh_te
    yerr = sqrt(xerr^2 + hii[indx].zt_log12oh_te_err^2)

    djs_plot, [0], [0], xrange=xrange, yrange=yrange, xtickname=replicate(' ',10), $
      charsize=charsize_6, charthick=postthick, xthick=postthick, ythick=postthick, $
      xsty=3, ysty=3, position=pos[*,0]
    oplot, x, y, ps=8
;   oploterror, 9.5*[1,1], -0.8*[1,1], djs_mean(xerr), djs_mean(yerr), ps=3, errthick=postthick, /nohat
;   oploterror, x, y, xerr, yerr, ps=8, errthick=postthick, /nohat
    djs_oplot, !x.crange, [0,0], thick=postthick

    legend, 'CL01 - Case A', /left, /bottom, box=0, charsize=charsize_6, charthick=postthick

; --------------------------------------------------    
; MacGaugh (1991) based on the Contini et al. (2002) decision tree 
; --------------------------------------------------    

    indx = where(hii.zstrong_12oh_M91_contini gt -900.0 and hii.zt_log12oh_te gt -900.0,nindx)

    x = hii[indx].zstrong_12oh_M91_contini
;   xerr = hii[indx].zstrong_12oh_M91_contini_err
    xerr = x*0.0

    y = hii[indx].zstrong_12oh_M91_contini - hii[indx].zt_log12oh_te
    yerr = sqrt(xerr^2 + hii[indx].zt_log12oh_te_err^2)

    djs_plot, [0], [0], xrange=xrange, yrange=yrange, charsize=charsize_6, $
      charthick=postthick, xthick=postthick, ythick=postthick, xsty=3, ysty=3, $
      position=pos[*,1], /noerase, xtickname=replicate(' ',10), ytickname=replicate(' ',10)
    oplot, x, y, ps=8
;   oploterror, x, y, xerr, yerr, ps=8, errthick=postthick, /nohat
    djs_oplot, !x.crange, [0,0], thick=postthick

    legend, 'M91 - r23', /left, /bottom, box=0, charsize=charsize_6, charthick=postthick

; --------------------------------------------------    
; Kewley & Dopita (2002) [N II]/[O II] calibration
; --------------------------------------------------    
 
    indx = where(hii.zstrong_12oh_kd02_niioii gt -900.0 and hii.zt_log12oh_te gt -900.0,nindx)

    x = hii[indx].zstrong_12oh_kd02_niioii
;   xerr = hii[indx].zstrong_12oh_kd02_niioii_err
    xerr = x*0.0

    y = hii[indx].zstrong_12oh_kd02_niioii-hii[indx].zt_log12oh_te
    yerr = sqrt(xerr^2 + hii[indx].zt_log12oh_te_err^2)

    djs_plot, [0], [0], xrange=xrange, yrange=yrange, charsize=charsize_6, $
      charthick=postthick, xthick=postthick, ythick=postthick, xsty=3, ysty=3, $
      position=pos[*,2], /noerase;, xtitle=xtitle
    oplot, x, y, ps=8
;   oploterror, x, y, xerr, yerr, ps=8, errthick=postthick, /nohat
    djs_oplot, !x.crange, [0,0], thick=postthick
 
    legend, 'KD02 - [N II]/[O II]', /left, /bottom, box=0, charsize=charsize_6, charthick=postthick

; --------------------------------------------------    
; Kewley & Dopita (2002) combined calibration
; --------------------------------------------------    

    indx = where(hii.zstrong_12oh_kd02_combined gt -900.0 and hii.zt_log12oh_te gt -900.0,nindx)

    x = hii[indx].zstrong_12oh_kd02_combined
;   xerr = hii[indx].zstrong_12oh_kd02_combined_err
    xerr = x*0.0

    y = hii[indx].zstrong_12oh_kd02_combined-hii[indx].zt_log12oh_te
    yerr = sqrt(xerr^2 + hii[indx].zt_log12oh_te_err^2)

    djs_plot, [0], [0], xrange=xrange, yrange=yrange, charsize=charsize_6, $
      charthick=postthick, xthick=postthick, ythick=postthick, xsty=3, ysty=3, $
      position=pos[*,3], /noerase, ytickname=replicate(' ',10);, xtitle=xtitle
    oplot, x, y, ps=8
;   oploterror, x, y, xerr, yerr, ps=8, errthick=postthick, /nohat
    djs_oplot, !x.crange, [0,0], thick=postthick

    legend, 'KD02 - Combined', /left, /bottom, box=0, charsize=charsize_6, charthick=postthick

; the titles

    xyouts, 0.05, 0.52, textoidl(ytitle), orientation=90, align=0.5, $
      charsize=charsize_6, charthick=postthick, /normal
    xyouts, 0.58, 0.025, textoidl(xtitle), align=0.5, $
      charsize=charsize_6, charthick=postthick, /normal

    im_openclose, postscript=postscript, /close    

;;; ------------------------------------------------------------
;;; [N II]/Ha vs Te
;;; ------------------------------------------------------------
;;
;;    psname = 'niiha_vs_te'
;;    im_openclose, pspath+psname, postscript=postscript
;;
;;    indx = where((hii.nii_6584_h_alpha gt -900.0) and (hii.zt_t4363 gt -900.0),nindx)
;;    
;;    x = hii[indx].nii_6584_h_alpha
;;    xerr = hii[indx].nii_6584_h_alpha_err
;;
;;    y = hii[indx].zt_t4363
;;    yerr = hii[indx].zt_t4363_err
;;
;;    xrange = [-2.6,-0.1]
;;    yrange = Terange
;;
;;    xtitle = 'log ([N II] \lambda6584 / H\alpha)'
;;    ytitle = 'T_{e} [K]'
;;
;;    zindicators_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
;;      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
;;      xstyle=3, ystyle=3, /right, /top, ytickname=Telabel
;;
;;; fit a linear function
;;
;;    sixlin, x, y, a, siga, b, sigb
;;    a = a[2] & b = b[2] & sigma_a_b = [siga[2],sigb[2]]
;;;   fitexy, x, y, a, b, x_sigma=xerr, y_sigma=yerr, sigma_a_b, chisq, q, tol=tol
;;    xaxis = findgen((max(x)-min(x))/0.001+1)*0.001+min(x)
;;    yfit = poly(xaxis,[a,b])
;;    djs_oplot, xaxis, yfit, line=0, thick=postthick, color='dark green'
;;
;;    splog, '[N II]/Ha vs Te coefficients: ', [a,b], sigma_a_b
;;
;;    legend, 'This Paper', line=0, thick=postthick, /right, /top, box=0, $
;;      charsize=1.5, charthick=postthick, color=djs_icolor('dark green')
;;    
;;    im_openclose, postscript=postscript, /close    
;;    
; ------------------------------------------------------------
; Te vs Te [empirical]
; ------------------------------------------------------------
;
;   psname = 'Te_direct_vs_te_niiha'
;   im_openclose, pspath+psname, postscript=postscript
;
;   indx = where((hii.zt_t4363 gt -900.0) and (hii.nii_6584_h_alpha gt -900.0),nindx)
;   
;   x = hii[indx].zt_t4363
;   xerr = hii[indx].zt_t4363_err
;
;   y = -7470.0*hii[indx].nii_6584_h_alpha + 2090
;   yerr = sqrt( (-7470.0*hii[indx].nii_6584_h_alpha_err)^2.0 + (64.0 * hii[indx].nii_6584_h_alpha)^2.0 + 101^2.0 )
;
;   xrange = Terange
;   yrange = Terange
;
;   xtitle = 'T_{e} [Direct, K]'
;   ytitle = 'T_{e} [Empirical, K]'
;
;   zindicators_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
;     xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
;     xstyle=3, ystyle=3, /right, /top, xtickname=Telabel, ytickname=Telabel
;   
;   im_openclose, postscript=postscript, /close    
;    
;;; ------------------------------------------------------------
;;; 12+log (O/H) [Te] vs ([N II]/Ha) / ([O III]/[O II])
;;; ------------------------------------------------------------
;;
;;    psname = '12oh_te_vs_niiha_oiiioii'
;;    im_openclose, pspath+psname, postscript=postscript
;;
;;    indx = where((hii.nii_6584_h_alpha gt -900.0) and (hii.zt_log12oh_te gt -900.0) and $
;;      (hii.oiii_5007_oii gt -900.0),nindx)
;;    
;;    x = hii[indx].zt_log12oh_te
;;    xerr = hii[indx].zt_log12oh_te_err
;;
;;    y = hii[indx].nii_6584_h_alpha - hii[indx].oiii_5007_oii
;;    yerr = sqrt(hii[indx].nii_6584_h_alpha_err^2.0 + hii[indx].oiii_5007_oii_err^2.0)
;;
;;    xrange = ohrange
;;    yrange = [-3.8,0.8]
;;
;;    xtitle = '12 + log (O/H) [T_{e}]'
;;    ytitle = 'log ([N II] \lambda6584/H\alpha)/([O III] \lambda5007/[O II] \lambda3727)'
;;    
;;    zindicators_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
;;      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
;;      /left, xstyle=3, ystyle=3, /right, /top, charsize=1.8
;;
;;    im_openclose, postscript=postscript, /close    
;;
;;; ------------------------------------------------------------
;;; 12+log (O/H) [Te] vs ([N II]/Ha) / ([O III]/Hb)
;;; ------------------------------------------------------------
;;    
;;    psname = '12oh_te_vs_niiha_oiiihb'
;;    im_openclose, pspath+psname, postscript=postscript
;;
;;    indx = where((hii.nii_6584_h_alpha gt -900.0) and (hii.zt_log12oh_te gt -900.0) and $
;;      (hii.oiii_5007_h_beta gt -900.0),nindx)
;;    
;;    x = hii[indx].zt_log12oh_te
;;    xerr = hii[indx].zt_log12oh_te_err
;;
;;    y = hii[indx].nii_6584_h_alpha - hii[indx].oiii_5007_h_beta
;;    yerr = sqrt(hii[indx].nii_6584_h_alpha_err^2.0 + hii[indx].oiii_5007_h_beta_err^2.0)
;;
;;    xrange = ohrange
;;    yrange = [-3.4,0.9]
;;
;;    xtitle = '12 + log (O/H) [T_{e}]'
;;    ytitle = 'log ([N II] \lambda6584/H\alpha)/([O III] \lambda5007/H\beta)'
;;    
;;    zindicators_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
;;      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
;;      /left, xstyle=3, ystyle=3, /right, /top, charsize=1.8
;;
;;    im_openclose, postscript=postscript, /close    
;;
;;; ------------------------------------------------------------
;;; [O III] / [O II] vs 12+log (O/H) Residuals
;;; ------------------------------------------------------------
;;
;;    psname = 'oiiioii_vs_12oh_residuals'
;;    im_openclose, pspath+psname, postscript=postscript
;;
;;    indx = where((hii.oiii_5007_oii gt -900.0) and (hii.zt_log12oh_te gt -900.0) and $
;;      (hii.nii_6584_h_alpha gt -900.0),nindx)
;;
;;    x = hii[indx].oiii_5007_oii
;;    xerr = hii[indx].oiii_5007_oii_err
;;
;;    oh_axis = hii[indx].nii_6584_h_alpha
;;    oh_model = a + b * oh_axis
;;    
;;    y = hii[indx].zt_log12oh_te - oh_model
;;    yerr = hii[indx].zt_log12oh_te_err
;;
;;    xrange = [-1.5,1.5]
;;    yrange = [-1.0,1.0]
;;
;;    xtitle = 'log ([O III] \lambda5007 / [O II] \lambda3727)'
;;    ytitle = '(12 + log (O/H) [T_{e}]) - (12 + log (O/H) [N II]/H\alpha)'
;;    
;;    zindicators_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
;;      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
;;      /left, xstyle=3, ystyle=3, /right, /top, charsize=1.8
;;    djs_oplot, !x.crange, [0,0], line=0, thick=postthick
;;    
;;    im_openclose, postscript=postscript, /close    
;;
;;; ------------------------------------------------------------
;;; [O III] / Hb vs 12+log (O/H) Residuals
;;; ------------------------------------------------------------
;;
;;    psname = 'oiiihb_vs_12oh_residuals'
;;    im_openclose, pspath+psname, postscript=postscript
;;
;;    indx = where((hii.oiii_5007_h_beta gt -900.0) and (hii.zt_log12oh_te gt -900.0) and $
;;      (hii.nii_6584_h_alpha gt -900.0),nindx)
;;
;;    x = hii[indx].oiii_5007_h_beta
;;    xerr = hii[indx].oiii_5007_h_beta_err
;;
;;    oh_axis = hii[indx].nii_6584_h_alpha
;;    oh_model = a + b * oh_axis
;;    
;;    y = hii[indx].zt_log12oh_te - oh_model
;;    yerr = hii[indx].zt_log12oh_te_err
;;
;;    xrange = [-1.5,1.2]
;;    yrange = [-1.0,1.0]
;;
;;    xtitle = 'log ([O III] \lambda5007 / H\beta)'
;;    ytitle = '(12 + log (O/H) [T_{e}]) - (12 + log (O/H) [N II]/H\alpha)'
;;    
;;    zindicators_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
;;      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
;;      /left, xstyle=3, ystyle=3, /right, /top, charsize=1.8
;;    djs_oplot, !x.crange, [0,0], line=0, thick=postthick
;;    
;;    im_openclose, postscript=postscript, /close    
;;
;;; ------------------------------------------------------------
;;; r23 vs 12+log (O/H) Residuals
;;; ------------------------------------------------------------
;;
;;    psname = 'r23_vs_12oh_residuals'
;;    im_openclose, pspath+psname, postscript=postscript
;;
;;    indx = where((hii.zstrong_r23 gt -900.0) and (hii.zt_log12oh_te gt -900.0) and $
;;      (hii.nii_6584_h_alpha gt -900.0),nindx)
;;
;;    x = alog10(hii[indx].zstrong_r23)
;;    xerr = hii[indx].zstrong_r23_err/alog10(hii[indx].zstrong_r23)/alog(10.0)
;;
;;    oh_axis = hii[indx].nii_6584_h_alpha
;;    oh_model = a + b * oh_axis
;;    
;;    y = hii[indx].zt_log12oh_te - oh_model
;;    yerr = hii[indx].zt_log12oh_te_err
;;
;;    xrange = [-0.5,1.2]
;;    yrange = [-1.0,1.0]
;;
;;    xtitle = 'log (R_{23})'
;;    ytitle = '(12 + log (O/H) [T_{e}]) - (12 + log (O/H) [N II]/H\alpha)'
;;    
;;    zindicators_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
;;      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
;;      /left, xstyle=3, ystyle=3, /right, /top, charsize=1.8
;;    djs_oplot, !x.crange, [0,0], line=0, thick=postthick
;;    
;;    im_openclose, postscript=postscript, /close    
;;
;;; ------------------------------------------------------------
;;; 12+log (O/H) [empirical] vs [N II]/Ha
;;; ------------------------------------------------------------
;;
;;    psname = '12oh_empirical_vs_niiha'
;;    im_openclose, pspath+psname, postscript=postscript
;;
;;    indx = where((hii.nii_6584_h_alpha gt -900.0) and (hii.zstrong_log12oh gt -900.0),nindx)
;;    
;;    x = hii[indx].zstrong_log12oh
;;    xerr = hii[indx].zstrong_log12oh_err
;;
;;    y = hii[indx].nii_6584_h_alpha
;;    yerr = hii[indx].nii_6584_h_alpha_err
;;
;;    xrange = ohrange
;;    yrange = [-2.6,-0.1]
;;
;;    xtitle = '12 + log (O/H) [empirical]'
;;    ytitle = 'log ([N II] \lambda6584 / H\alpha)'
;;    
;;    zindicators_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
;;      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
;;      /left, xstyle=3, ystyle=3, /right, /top
;;
;;    im_openclose, postscript=postscript, /close    
;;
;;; ------------------------------------------------------------
;;; 12+log (O/H) [empirical] vs [N II]/Ha - Te sample only
;;; ------------------------------------------------------------
;;
;;    psname = '12oh_empirical_vs_niiha_te_subsample'
;;    im_openclose, pspath+psname, postscript=postscript
;;
;;    indx = where((hii.nii_6584_h_alpha gt -900.0) and (hii.zstrong_log12oh gt -900.0) and $
;;      (hii.zt_t_oiii gt -900.0),nindx)
;;    
;;    x = hii[indx].zstrong_log12oh
;;    xerr = hii[indx].zstrong_log12oh_err
;;
;;    y = hii[indx].nii_6584_h_alpha
;;    yerr = hii[indx].nii_6584_h_alpha_err
;;
;;    xrange = ohrange
;;    yrange = [-2.6,-0.1]
;;
;;    xtitle = '12 + log (O/H) [empirical]'
;;    ytitle = 'log ([N II] \lambda6584 / H\alpha)'
;;    
;;    zindicators_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
;;      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
;;      /left, xstyle=3, ystyle=3, /right, /top
;;
;;    legend, textoidl('T_{e} Subsample'), /right, /bottom, box=0, $
;;      charsize=1.5, charthick=postthick
;;
;;    im_openclose, postscript=postscript, /close    
;;
;;; ------------------------------------------------------------
;;; 12+log (O/H) [niiha] vs [N II]/[O II]
;;; ------------------------------------------------------------
;;
;;    psname = '12oh_niiha_vs_niioii'
;;    im_openclose, pspath+psname, postscript=postscript
;;
;;    indx = where((hii.nii_6584_oii gt -900.0) and (hii.zstrong_12oh_n2 gt -900.0),nindx)
;;    
;;    x = hii[indx].zstrong_12oh_n2
;;    xerr = hii[indx].zstrong_12oh_n2_err
;;
;;    y = hii[indx].nii_6584_oii
;;    yerr = hii[indx].nii_6584_oii_err
;;
;;    xrange = ohrange
;;    yrange = [-1.8,0.6]
;;
;;    xtitle = '12 + log (O/H) [N II]/H\alpha'
;;    ytitle = 'log ([N II] \lambda6584 / [O II] \lambda3727)'
;;    
;;    zindicators_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
;;      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
;;      /left, xstyle=3, ystyle=3, /right, /top
;;       
;;    im_openclose, postscript=postscript, /close    
;;
;;; ------------------------------------------------------------
;;; 12+log (O/H) [empirical] vs [O III]/Hb
;;; ------------------------------------------------------------
;;
;;    psname = '12oh_empirical_vs_oiiihb'
;;    im_openclose, pspath+psname, postscript=postscript
;;
;;    indx = where((hii.oiii_5007_h_beta gt -900.0) and (hii.zstrong_log12oh gt -900.0),nindx)
;;    
;;    x = hii[indx].zstrong_log12oh
;;    xerr = hii[indx].zstrong_log12oh_err
;;
;;    y = hii[indx].oiii_5007_h_beta
;;    yerr = hii[indx].oiii_5007_h_beta_err
;;
;;    xrange = ohrange
;;    yrange = [-1.5,1.2]
;;
;;    xtitle = '12 + log (O/H) [empirical]'
;;    ytitle = 'log ([O III] \lambda5007 / H\beta)'
;;
;;    zindicators_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
;;      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
;;      xstyle=3, ystyle=3, /right, /top
;;
;;    im_openclose, postscript=postscript, /close    
;;
;;; ------------------------------------------------------------
;;; 12+log (O/H) [niiha] vs [O III]/Hb
;;; ------------------------------------------------------------
;;
;;    psname = '12oh_niiha_vs_oiiihb'
;;    im_openclose, pspath+psname, postscript=postscript
;;
;;    indx = where((hii.oiii_5007_h_beta gt -900.0) and (hii.zstrong_12oh_n2 gt -900.0),nindx)
;;    
;;    x = hii[indx].zstrong_12oh_n2
;;    xerr = hii[indx].zstrong_12oh_n2_err
;;
;;    y = hii[indx].oiii_5007_h_beta
;;    yerr = hii[indx].oiii_5007_h_beta_err
;;
;;    xrange = ohrange
;;    yrange = [-1.5,1.2]
;;
;;    xtitle = '12 + log (O/H) [N II]/H\alpha'
;;    ytitle = 'log ([O III] \lambda5007 / H\beta)'
;;
;;    zindicators_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
;;      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
;;      xstyle=3, ystyle=3, /right, /top
;;
;;    im_openclose, postscript=postscript, /close    
;;    
;;; ------------------------------------------------------------
;;; 12+log (O/H) [niiha] vs r23
;;; ------------------------------------------------------------
;;
;;    psname = '12oh_niiha_vs_r23'
;;    im_openclose, pspath+psname, postscript=postscript
;;
;;    indx = where((hii.zstrong_r23 gt -900.0) and (hii.zstrong_12oh_n2 gt -900.0),nindx)
;;    
;;    x = hii[indx].zstrong_12oh_n2
;;    xerr = hii[indx].zstrong_12oh_n2_err
;;
;;    y = alog10(hii[indx].zstrong_r23)
;;    yerr = hii[indx].zstrong_r23_err/alog10(hii[indx].zstrong_r23)/alog(10.0)
;;
;;    xrange = ohrange
;;    yrange = [-0.5,1.2]
;;
;;    xtitle = '12 + log (O/H) [N II]/H\alpha'
;;    ytitle = 'log (R_{23})'
;;    
;;    zindicators_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
;;      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
;;      /left, xstyle=3, ystyle=3, /right, /top
;;
;;    im_openclose, postscript=postscript, /close    
;;
;;; ------------------------------------------------------------
;;; r23 vs 12+log (O/H)
;;; ------------------------------------------------------------
;;
;;    psname = 'r23_vs_12oh'
;;    im_openclose, pspath+psname, postscript=postscript
;;
;;    indx = where((hii.zstrong_r23 gt -900.0) and (hii.zstrong_12oh_n2 gt -900.0),nindx)
;;     
;;    x = alog10(hii[indx].zstrong_r23)
;;    xerr = hii[indx].zstrong_r23_err/alog10(hii[indx].zstrong_r23)/alog(10.0)
;;
;;    y = hii[indx].zstrong_12oh_n2
;;    yerr = hii[indx].zstrong_12oh_n2_err
;;
;;; -----    
;;    indx1 = where((hii.zstrong_r23 gt -900.0) and (hii.zt_log12oh_te gt -900.0),nindx1)
;;     
;;    x1 = hii[indx1].zstrong_r23
;;    x1err = hii[indx1].zstrong_r23_err
;;
;;    y1 = hii[indx1].zt_log12oh_te
;;    y1err = hii[indx1].zt_log12oh_te_err
;;; -----    
;;
;;    xrange = [-0.5,1.2]
;;    yrange = ohrange
;;
;;    xtitle = 'log (R_{23})'
;;    ytitle = '12 + log (O/H)'
;;    
;;    zindicators_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
;;      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
;;      /left, xstyle=3, ystyle=3, /right, /top
;;    plotsym, 8, 0.8, /fill
;;;   oploterror, x1, y1, x1err, y1err, ps=8, /nohat, thick=postthick, $
;;;     errthick=postthick, color=djs_icolor('red'), errcolor=djs_icolor('red')
;;    oplot, x1, y1, ps=8, thick=postthick, color=djs_icolor('red')
;;
;;    im_legend, textoidl(['[N II]/H\alpha','T_{e}']), /left, /bottom, box=0, $
;;      charsize=2.0, charthick=postthick, color=djs_icolor(['purple','red']), $
;;      psym=[108,106], fill=[1,1], symsize=1.3
;;
;;    im_openclose, postscript=postscript, /close    
;;
;;; ------------------------------------------------------------
;;; [O III]/[O II] vs 12+log (O/H)
;;; ------------------------------------------------------------
;;
;;    psname = 'oiiioii_vs_12oh'
;;    im_openclose, pspath+psname, postscript=postscript
;;
;;    indx = where((hii.oiii_5007_oii gt -900.0) and (hii.zstrong_12oh_n2 gt -900.0),nindx)
;;    
;;    x = hii[indx].oiii_5007_oii
;;    xerr = hii[indx].oiii_5007_oii_err
;;
;;    y = hii[indx].zstrong_12oh_n2
;;    yerr = hii[indx].zstrong_12oh_n2_err
;;
;;; -----    
;;    indx1 = where((hii.oiii_5007_oii gt -900.0) and (hii.zt_log12oh_te gt -900.0),nindx1)
;;     
;;    x1 = hii[indx1].oiii_5007_oii
;;    x1err = hii[indx1].oiii_5007_oii_err
;;
;;    y1 = hii[indx1].zt_log12oh_te
;;    y1err = hii[indx1].zt_log12oh_te_err
;;; -----    
;;
;;    xrange = [-1.5,1.5]
;;    yrange = ohrange
;;
;;    xtitle = 'log ([O III] \lambda5007 / [O II] \lambda3727)'
;;    ytitle = '12 + log (O/H)'
;;    
;;    zindicators_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
;;      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
;;      /left, xstyle=3, ystyle=3, /right, /top
;;    plotsym, 8, 0.8, /fill
;;;   oploterror, x1, y1, x1err, y1err, ps=8, /nohat, thick=postthick, $
;;;     errthick=postthick, color=djs_icolor('red'), errcolor=djs_icolor('red')
;;    oplot, x1, y1, ps=8, thick=postthick, color=djs_icolor('red')
;;
;;    im_legend, textoidl(['[N II]/H\alpha','T_{e}']), /right, /top, box=0, $
;;      charsize=2.0, charthick=postthick, color=djs_icolor(['purple','red']), $
;;      psym=[108,106], fill=[1,1], symsize=1.3
;;
;;    im_openclose, postscript=postscript, /close    
;;
;;; ------------------------------------------------------------
;;; 12 + log (O/H) [Empirical] vs 12 + log (O/H) [Te]
;;; ------------------------------------------------------------
;;
;;    psname = '12oh_empirical_vs_12oh_te'
;;    im_openclose, pspath+psname, postscript=postscript
;;
;;    indx = where((hii.zt_log12oh_te gt -900.0) and (hii.zstrong_log12oh gt -900.0),nindx)
;;    
;;    x = hii[indx].zstrong_log12oh
;;    xerr = hii[indx].zstrong_log12oh_err
;;
;;    y = hii[indx].zt_log12oh_te
;;    yerr = hii[indx].zt_log12oh_te_err
;;
;;    xrange = ohrange
;;    yrange = ohrange
;;
;;    xtitle = '12 + log (O/H) [Empirical]'
;;    ytitle = '12 + log (O/H) [T_{e}]'
;;
;;    zindicators_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
;;      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
;;      xstyle=3, ystyle=3, /right, /top
;;    djs_oplot, !x.crange, !y.crange, thick=postthick
;;
;;    im_openclose, postscript=postscript, /close    
;;
;;; ------------------------------------------------------------
;;; 12 + log (O/H) [Te] vs 12 + log (O/H) [P01] - Plot I
;;; ------------------------------------------------------------
;;
;;    psname = '12oh_te_vs_12oh_residuals'
;;    im_openclose, pspath+psname, postscript=postscript
;;
;;    lower = where((hii.zstrong_12oh_P01_lower gt -900.0) and (hii.zt_log12oh_te gt -900.0) and $
;;      (hii.zt_log12oh_te le 7.95),nlower)
;;    upper = where((hii.zstrong_12oh_P01_upper gt -900.0) and (hii.zt_log12oh_te gt -900.0) and $
;;      (hii.zt_log12oh_te ge 8.2),nupper)
;;
;;    indx = cmset_op(lower,'OR',upper)
;;
;;; -----
;;    Plo = where(hii[lower].zstrong_P lt 0.4,comp=Phi)
;;
;;;   xlower = hii[lower].zt_log12oh_te
;;;   xlowererr = hii[lower].zt_log12oh_te_err
;;
;;    xlower_Plo = hii[lower[Plo]].zt_log12oh_te
;;    xlowererr_Plo = hii[lower[Plo]].zt_log12oh_te_err
;;
;;    xlower_Phi = hii[lower[Phi]].zt_log12oh_te
;;    xlowererr_Phi = hii[lower[Phi]].zt_log12oh_te_err
;;
;;;   ylower = hii[lower].zstrong_12oh_P01_lower - xlower
;;;   ylowererr = sqrt(hii[lower].zstrong_12oh_P01_lower_err^2 + xlowererr^2)
;;
;;    ylower_Plo = hii[lower[Plo]].zstrong_12oh_P01_lower - xlower_Plo
;;    ylowererr_Plo = sqrt(hii[lower[Plo]].zstrong_12oh_P01_lower_err^2 + xlowererr_Plo^2)
;;
;;    ylower_Phi = hii[lower[Phi]].zstrong_12oh_P01_lower - xlower_Phi
;;    ylowererr_Phi = sqrt(hii[lower[Phi]].zstrong_12oh_P01_lower_err^2 + xlowererr_Phi^2)
;;; -----
;;
;;    Plo = where(hii[upper].zstrong_P lt 0.4,comp=Phi)
;;
;;;   xupper = hii[upper].zt_log12oh_te
;;;   xuppererr = hii[upper].zt_log12oh_te_err
;;
;;    xupper_Plo = hii[upper[Plo]].zt_log12oh_te
;;    xuppererr_Plo = hii[upper[Plo]].zt_log12oh_te_err
;;
;;    xupper_Phi = hii[upper[Phi]].zt_log12oh_te
;;    xuppererr_Phi = hii[upper[Phi]].zt_log12oh_te_err
;;
;;;   yupper = hii[upper].zstrong_12oh_P01_upper - xupper
;;;   yuppererr = sqrt(hii[upper].zstrong_12oh_P01_upper_err^2 + xuppererr^2)
;;
;;    yupper_Plo = hii[upper[Plo]].zstrong_12oh_P01_upper - xupper_Plo
;;    yuppererr_Plo = sqrt(hii[upper[Plo]].zstrong_12oh_P01_upper_err^2 + xuppererr_Plo^2)
;;
;;    yupper_Phi = hii[upper[Phi]].zstrong_12oh_P01_upper - xupper_Phi
;;    yuppererr_Phi = sqrt(hii[upper[Phi]].zstrong_12oh_P01_upper_err^2 + xuppererr_Phi^2)
;;; -----
;;
;;    xrange = ohrange
;;    yrange = [-0.75,0.75]
;;
;;    xtitle = '12 + log (O/H) [T_{e}]'
;;    ytitle = 'log (O/H) [P01] - log (O/H) T_{e}'
;;
;;    djs_plot, [0], [0], /nodata, xtitle=xtitle, ytitle=ytitle, $
;;      xrange=xrange, yrange=yrange, xstyle=3, ystyle=3, $
;;      xthick=postthick, ythick=postthick, charsize=2.0, charthick=postthick
;;    djs_oplot, !x.crange, [0,0], thick=postthick
;;
;;    plotsym, 0, 1, /fill
;;    djs_oplot, [xlower_Phi,xupper_Phi], [ylower_Phi,yupper_Phi], ps=8, color='red'
;;
;;    plotsym, 8, 1, /fill
;;    djs_oplot, [xlower_Plo,xupper_Plo], [ylower_Plo,yupper_Plo], ps=8, color='green'
;;    
;;    im_legend, ['P > 0.4','P < 0.4'], /right, /top, box=0, charsize=2.0, $
;;      charthick=postthick, psym=[108,106], fill=[1,1], symsize=1.4, $
;;      color=djs_icolor(['red','green'])
;;    
;;    im_openclose, postscript=postscript, /close    

;;; ------------------------------------------------------------
;;; Various Abundances vs Te-based abundances - Part II
;;; ------------------------------------------------------------
;;
;;    psname = 'ZTe_vs_Zstrong2'
;;    im_openclose, pspath+psname, postscript=postscript
;;
;;    plotsym, 8, 0.5, /fill
;;
;;    pagemaker, nx=2, ny=4, position=pos, /normal, xspace=0.0, $
;;      xmargin=[1.6,0.1], ymargin=[0.5,1.2], yspace=0.0
;;
;;    ytitle = '12 + log (O/H) [Strong-Line]'
;;    xtitle = '12 + log (O/H) [T_{e}]'
;;
;;    xrange = [7.0,9.9]
;;    yrange = [7.0,9.9]
;;    charsize_6 = 1.2
;;    charsize_6 = 1.4
;;
;;; --------------------------------------------------    
;;; Charlot & Longhetti (2001) - Case A
;;; --------------------------------------------------    
;;
;;    indx = where(hii.zstrong_12oh_CL01A gt -900.0 and hii.zt_log12oh_te gt -900.0,nindx)
;;
;;    x = hii[indx].zt_log12oh_te
;;    xerr = hii[indx].zt_log12oh_te_err
;;
;;    y = hii[indx].zstrong_12oh_CL01A
;;    yerr = hii[indx].zstrong_12oh_CL01A_err
;;
;;    djs_plot, [0], [0], xrange=xrange, yrange=yrange, xtickname=replicate(' ',10), $
;;      charsize=charsize_6, charthick=postthick, xthick=postthick, ythick=postthick, $
;;      xsty=3, ysty=3, position=pos[*,0]
;;    oplot, x, y, ps=8
;;;   oploterror, x, y, xerr, yerr, ps=8, errthick=postthick, /nohat
;;    djs_oplot, !x.crange, !y.crange, thick=postthick
;;
;;    legend, 'CL01 - Case A', /right, /bottom, box=0, charsize=charsize_6, charthick=postthick
;;
;;; --------------------------------------------------    
;;; Charlot & Longhetti (2001) - Case F
;;; --------------------------------------------------    
;;
;;    indx = where(hii.zstrong_12oh_CL01F gt -900.0 and hii.zt_log12oh_te gt -900.0,nindx)
;;
;;    x = hii[indx].zt_log12oh_te
;;    xerr = hii[indx].zt_log12oh_te_err
;;
;;    y = hii[indx].zstrong_12oh_CL01F
;;    yerr = hii[indx].zstrong_12oh_CL01F_err
;;
;;    djs_plot, [0], [0], xrange=xrange, yrange=yrange, charsize=charsize_6, charthick=postthick, $
;;      xthick=postthick, ythick=postthick, xsty=3, ysty=3, xtickname=replicate(' ',10), $
;;      ytickname=replicate(' ',10), position=pos[*,1], /noerase
;;    oplot, x, y, ps=8
;;;   oploterror, x, y, xerr, yerr, ps=8, errthick=postthick, /nohat
;;    djs_oplot, !x.crange, !y.crange, thick=postthick
;;
;;    legend, 'CL01 - Case F', /right, /bottom, box=0, charsize=charsize_6, charthick=postthick
;;
;;; --------------------------------------------------    
;;; Zaritsky et al. (1994)
;;; --------------------------------------------------    
;;
;;    indx = where(hii.zstrong_12oh_ZKH94 gt -900.0 and hii.zt_log12oh_te gt -900.0,nindx)
;;
;;    x = hii[indx].zt_log12oh_te
;;    xerr = hii[indx].zt_log12oh_te_err
;;
;;    y = hii[indx].zstrong_12oh_ZKH94
;;;   yerr = hii[indx].zstrong_12oh_ZKH94_err
;;    yerr = y*0.0+0.1
;;
;;    djs_plot, [0], [0], xrange=xrange, yrange=yrange, charsize=charsize_6, $
;;      charthick=postthick, xthick=postthick, ythick=postthick, xsty=3, ysty=3, $
;;      position=pos[*,2], /noerase, xtickname=replicate(' ',10)
;;    oplot, x, y, ps=8
;;;   oploterror, x, y, xerr, yerr, ps=8, errthick=postthick, /nohat
;;    djs_oplot, !x.crange, !y.crange, thick=postthick
;;    djs_oplot, 8.4*[1,1], !y.crange, line=2, thick=postthick
;;
;;    legend, 'ZKH94 - r23', /right, /bottom, box=0, charsize=charsize_6, charthick=postthick
;;
;;; --------------------------------------------------    
;;; Pilyugin (2001) & Pilyugin (2003) based on the Melbourne & Salzer
;;; (2002) decision tree
;;; --------------------------------------------------    
;;
;;    indx = where(hii.zstrong_12oh_P01_melbourne gt -900.0 and hii.zt_log12oh_te gt -900.0,nindx)
;;
;;    x = hii[indx].zt_log12oh_te
;;    xerr = hii[indx].zt_log12oh_te_err
;;
;;    y = hii[indx].zstrong_12oh_P01_melbourne
;;    yerr = hii[indx].zstrong_12oh_P01_melbourne_err
;;
;;    djs_plot, [0], [0], xrange=xrange, yrange=yrange, charsize=charsize_6, $
;;      charthick=postthick, xthick=postthick, ythick=postthick, xsty=3, ysty=3, $
;;      position=pos[*,3], /noerase, ytickname=replicate(' ',10), xtickname=replicate(' ',10)
;;    oplot, x, y, ps=8
;;;   oploterror, x, y, xerr, yerr, ps=8, errthick=postthick, /nohat
;;    djs_oplot, !x.crange, !y.crange, thick=postthick
;;
;;    legend, 'MS02 - r23', /right, /bottom, box=0, charsize=charsize_6, charthick=postthick
;;
;;; --------------------------------------------------    
;;; MacGaugh (1991) based on [O III]/Hb
;;; --------------------------------------------------    
;;
;;    indx = where(hii.zstrong_12oh_M91_o32 gt -900.0 and hii.zt_log12oh_te gt -900.0,nindx)
;;
;;    x = hii[indx].zt_log12oh_te
;;    xerr = hii[indx].zt_log12oh_te_err
;;
;;    y = hii[indx].zstrong_12oh_M91_o32
;;;   yerr = hii[indx].zstrong_12oh_M91_o32_err
;;    yerr = y*0.0+0.1
;;
;;    djs_plot, [0], [0], xrange=xrange, yrange=yrange, charsize=charsize_6, $
;;      charthick=postthick, xthick=postthick, ythick=postthick, xsty=3, ysty=3, $
;;      position=pos[*,4], /noerase, xtickname=replicate(' ',10)
;;    oplot, x, y, ps=8
;;;   oploterror, x, y, xerr, yerr, ps=8, errthick=postthick, /nohat
;;    djs_oplot, !x.crange, !y.crange, thick=postthick
;;
;;    legend, 'M91 - r23', /right, /bottom, box=0, charsize=charsize_6, charthick=postthick
;;
;;;; --------------------------------------------------    
;;;; MacGaugh (1991) based on the Contini et al. (2002) decision tree 
;;;; --------------------------------------------------    
;;;;
;;;;   indx = where(hii.zstrong_12oh_M91_contini gt -900.0 and hii.zt_log12oh_te gt -900.0,nindx)
;;;;
;;;;   x = hii[indx].zt_log12oh_te
;;;;   xerr = hii[indx].zt_log12oh_te_err
;;;;
;;;;   y = hii[indx].zstrong_12oh_M91_contini
;;;;   yerr = hii[indx].zstrong_12oh_M91_contini_err
;;;;   yerr = y*0.0+0.1
;;;;
;;;;   djs_plot, [0], [0], xrange=xrange, yrange=yrange, charsize=charsize_6, $
;;;;     charthick=postthick, xthick=postthick, ythick=postthick, xsty=3, ysty=3, $
;;;;     position=pos[*,4], /noerase, xtickname=replicate(' ',10)
;;;;   oplot, x, y, ps=8
;;;;   oploterror, x, y, xerr, yerr, ps=8, errthick=postthick, /nohat
;;;;   djs_oplot, !x.crange, !y.crange, thick=postthick
;;;;
;;;;   legend, 'C02 - r23', /right, /bottom, box=0, charsize=charsize_6, charthick=postthick
;;
;;; --------------------------------------------------    
;;; Kewley & Dopita (2002) [N II]/[O II] calibration
;;; --------------------------------------------------    
;; 
;;    indx = where(hii.zstrong_12oh_kd02_niioii gt -900.0 and hii.zt_log12oh_te gt -900.0,nindx)
;;
;;    x = hii[indx].zt_log12oh_te
;;    xerr = hii[indx].zt_log12oh_te_err
;;
;;    y = hii[indx].zstrong_12oh_kd02_niioii
;;;   yerr = hii[indx].zstrong_12oh_kd02_niioii_err
;;    yerr = y*0.0+0.1
;;
;;    djs_plot, [0], [0], xrange=xrange, yrange=yrange, $
;;      charsize=charsize_6, charthick=postthick, xthick=postthick, ythick=postthick, xsty=3, ysty=3, $
;;      position=pos[*,5], /noerase, ytickname=replicate(' ',10), xtickname=replicate(' ',10)
;;    oplot, x, y, ps=8
;;;   oploterror, x, y, xerr, yerr, ps=8, errthick=postthick, /nohat
;;    djs_oplot, !x.crange, !y.crange, thick=postthick
;;;   djs_oplot, !x.crange, 8.6*[1,1], line=2, thick=postthick
;; 
;;    legend, 'KD02 - [N II]/[O II]', /right, /bottom, box=0, charsize=charsize_6, charthick=postthick
;;
;;; --------------------------------------------------    
;;; P01
;;; --------------------------------------------------    
;;
;;    indx = where(hii.zstrong_12oh_P01 gt -900.0 and hii.zt_log12oh_te gt -900.0,nindx)
;;
;;    x = hii[indx].zt_log12oh_te
;;    xerr = hii[indx].zt_log12oh_te_err
;;
;;    y = hii[indx].zstrong_12oh_P01
;;;   yerr = hii[indx].zstrong_12oh_P01_err
;;    yerr = y*0.0+0.1
;;
;;    djs_plot, [0], [0], xrange=xrange, yrange=yrange, xtitle=xtitle, charsize=charsize_6, $
;;      charthick=postthick, xthick=postthick, ythick=postthick, xsty=3, ysty=3, $
;;      position=pos[*,6], /noerase
;;    oplot, x, y, ps=8
;;;   oploterror, x, y, xerr, yerr, ps=8, errthick=postthick, /nohat
;;    djs_oplot, !x.crange, !y.crange, thick=postthick
;;
;;    legend, 'P01', /right, /bottom, box=0, charsize=charsize_6, charthick=postthick
;;
;;; --------------------------------------------------    
;;; KK04 - r23
;;; --------------------------------------------------    
;;
;;    indx = where(hii.zstrong_12oh_KK04_o32 gt -900.0 and hii.zt_log12oh_te gt -900.0,nindx)
;;
;;    x = hii[indx].zt_log12oh_te
;;    xerr = hii[indx].zt_log12oh_te_err
;;
;;    y = hii[indx].zstrong_12oh_KK04_o32
;;;   yerr = hii[indx].zstrong_12oh_KK04_o32_err
;;    yerr = y*0.0+0.1
;;
;;    djs_plot, [0], [0], xrange=xrange, yrange=yrange, xtitle=xtitle, charsize=charsize_6, $
;;      charthick=postthick, xthick=postthick, ythick=postthick, xsty=3, ysty=3, $
;;      position=pos[*,7], /noerase, ytickname=replicate(' ',10)
;;    oplot, x, y, ps=8
;;;   oploterror, x, y, xerr, yerr, ps=8, errthick=postthick, /nohat
;;    djs_oplot, !x.crange, !y.crange, thick=postthick
;;
;;    legend, 'KK04 - r23', /right, /bottom, box=0, charsize=charsize_6, charthick=postthick
;;
;;; the title
;;
;;    xyouts, 0.06, 0.52, ytitle, orientation=90, align=0.5, charsize=charsize_6, $
;;      charthick=postthick, /normal
;;
;;    im_openclose, postscript=postscript, /close    
;;
;;; ------------------------------------------------------------
;;; log U vs [O III]/[O II]
;;; ------------------------------------------------------------
;;
;;    psname = 'logU_vs_oiii_oii'
;;    im_openclose, pspath+psname, postscript=postscript
;;
;;    plot_kewley_grids, plotnumber=18, postscript=postscript
;;    
;;    im_openclose, postscript=postscript, /close    
;;    
;;; ------------------------------------------------------------
;;; log Z vs [O III]/[O II]
;;; ------------------------------------------------------------
;;
;;    psname = 'logzstrong_vs_oiii_oii'
;;    im_openclose, pspath+psname, postscript=postscript
;;
;;    plot_kewley_grids, plotnumber=19, label=5, postscript=postscript
;;    
;;    im_openclose, postscript=postscript, /close    
;

;;; ------------------------------------------------------------
;;; 12 + log (O/H) [Empirical] vs 12 + log (O/H) [niiha]
;;; ------------------------------------------------------------
;;
;;    psname = '12oh_empirical_vs_12oh_niiha'
;;    im_openclose, pspath+psname, postscript=postscript
;;
;;    indx = where((hii.zstrong_log12oh gt -900.0) and (hii.zstrong_12oh_n2 gt -900.0),nindx)
;;    
;;    x = hii[indx].zstrong_log12oh
;;    xerr = hii[indx].zstrong_log12oh_err
;;
;;    y = hii[indx].zstrong_12oh_n2
;;    yerr = hii[indx].zstrong_12oh_n2_err
;;
;;    xrange = ohrange
;;    yrange = ohrange
;;
;;    xtitle = '12 + log (O/H) [Empirical]'
;;    ytitle = '12 + log (O/H) [N II]/H\alpha'
;;
;;    zindicators_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
;;      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
;;      xstyle=3, ystyle=3, /right, /top
;;    djs_oplot, !x.crange, !y.crange, thick=postthick
;;    
;;    im_openclose, postscript=postscript, /close    
;;
;;; ------------------------------------------------------------
;;; 12 + log (O/H) [niiha] vs 12 + log (O/H) [P01]
;;; ------------------------------------------------------------
;;
;;    psname = '12oh_niiha_vs_12oh_P01'
;;    im_openclose, pspath+psname, postscript=postscript
;;
;;    indx = where((hii.zstrong_12oh_P01 gt -900.0) and (hii.zstrong_12oh_n2 gt -900.0),nindx)
;;    
;;    x = hii[indx].zstrong_12oh_P01
;;;   xerr = hii[indx].zstrong_12oh_P01_err
;;    xerr = x*0.0+0.1
;;
;;    y = hii[indx].zstrong_12oh_n2
;;    yerr = hii[indx].zstrong_12oh_n2_err
;;
;;    xrange = ohrange
;;    yrange = ohrange
;;
;;    xtitle = '12 + log (O/H) [N II]/H\alpha'
;;    ytitle = '12 + log (O/H) [P01]'
;;
;;    zindicators_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
;;      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
;;      xstyle=3, ystyle=3, /right, /top
;;    djs_oplot, !x.crange, !y.crange, thick=postthick
;;    
;;    im_openclose, postscript=postscript, /close    
;;

; ------------------------------------------------------------
; 12+log(O/H) [oiiinii+niiha] vs 12+log(O/H) [Strong]
; ------------------------------------------------------------

    psname = 'hii_12oh_oiiinii_niiha_vs_12oh_strong'
    im_openclose, pspath+psname, postscript=postscript

    pagemaker, nx=1, ny=3, position=pos, /normal, xspace=0.0, $
      xmargin=[1.7,0.6], ymargin=[0.2,1.1], yspace=0.0

    ytitle = 'log (O/H) [Empirical] - log (O/H) T_{e}'
    xtitle = '12 + log (O/H) ([O III]/H\beta)/([N II]/H\alpha)'

    xrange = [7.2,9.8]
    yrange = [7.2,9.8]
    charsize_6 = 1.3

; --------------------------------------------------    
; KD02-Combined
; --------------------------------------------------    

    indx = where((hii.zstrong_12oh_oiiinii_pettini gt -900) and (hii.zstrong_12oh_kd02_niioii gt -900.0),nindx)

    x = hii[indx].zstrong_12oh_oiiinii_pettini
    xerr = hii[indx].zstrong_12oh_oiiinii_pettini_err

    y = hii[indx].zstrong_12oh_kd02_niioii
    yerr = hii[indx].zstrong_12oh_kd02_niioii_err

    ytitle = '12 + log (O/H) [N II]/[O II]'

    stats = im_stats(x-y,/verbose,/no_head)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)
    
    zindicators_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=3, ystyle=3, /right, /top, position=pos[*,0], charsize=charsize_6, xtickname=replicate(' ',10)
    djs_oplot, !x.crange, !y.crange, thick=postthick

    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_6, charthick=postthick
    legend, '(a)', /left, /top, box=0, charsize=charsize_6, charthick=postthick
    
; --------------------------------------------------    
; M91
; --------------------------------------------------    

    indx = where((hii.zstrong_12oh_oiiinii_pettini gt -900) and (hii.zstrong_12oh_m91_upper gt -900.0) and $
      (hii.zstrong_12oh_m91_lower gt -900.0),nindx)

    x = hii[indx].zstrong_12oh_oiiinii_pettini
    xerr = hii[indx].zstrong_12oh_oiiinii_pettini_err
    
    y = x*0.0
    yerr = xerr*0.0
    
    for idiv = 0L, ndiv-1L do begin

       up = where((x gt div[idiv]),nup,comp=lo,ncomp=nlo)
       if (nup ne 0L) then begin
          y[up] = hii[indx[up]].zstrong_12oh_m91_upper
          yerr[up] = hii[indx[up]].zstrong_12oh_m91_upper_err
       endif

       if (nlo ne 0L) then begin
          y[lo] = hii[indx[lo]].zstrong_12oh_m91_lower
          yerr[lo] = hii[indx[lo]].zstrong_12oh_m91_lower_err
       endif

       resid[idiv] = stddev(x-y)
;      resid[idiv] = djsig(x-y,sigrej=3.0)
;      plot, x, x-y, ps=4, xrange=xrange, yrange=[-1,1]
;      print, div[idiv], resid[idiv]
;      cc = get_kbrd(1)
       
    endfor

;   mindiv = find_nminima(resid,div,nfind=1,width=0.2)
    minresid = min(resid,minindx)
    mindiv = div[minindx]
    splog, 'oiiinii+niiha vs M91: '+string(mindiv,format='(F4.2)')

    up = where((hii[indx].zstrong_12oh_oiiinii_pettini gt mindiv),nup)
    if (nup ne 0L) then begin
       x = hii[indx[up]].zstrong_12oh_oiiinii_pettini
       xerr = hii[indx[up]].zstrong_12oh_oiiinii_pettini_err

       y = hii[indx[up]].zstrong_12oh_m91_upper
       yerr = hii[indx[up]].zstrong_12oh_m91_upper_err
    endif

    lo = where((hii[indx].zstrong_12oh_oiiinii_pettini lt mindiv),nlo)
    if (nlo ne 0L) then begin
       x = [x,hii[indx[lo]].zstrong_12oh_oiiinii_pettini]
       xerr = [xerr,hii[indx[lo]].zstrong_12oh_oiiinii_pettini_err]

       y = [y,hii[indx[lo]].zstrong_12oh_m91_lower]
       yerr = [yerr,hii[indx[lo]].zstrong_12oh_m91_lower_err]
    endif

    ytitle = '12 + log (O/H) M91'

    stats = im_stats(x-y,/verbose,/no_head)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)
    
    zindicators_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle='', ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, $
      xstyle=3, ystyle=3, /right, /top, position=pos[*,1], /noerase, charsize=charsize_6, $
      xtickname=replicate(' ',10)
    djs_oplot, !x.crange, !y.crange, thick=postthick

    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_6, charthick=postthick
    legend, '(b)', /left, /top, box=0, charsize=charsize_6, charthick=postthick
    
; --------------------------------------------------    
; KD02-Combined
; --------------------------------------------------    

    indx = where((hii.zstrong_12oh_oiiinii_pettini gt -900.0) and (hii.zstrong_12oh_kd02_combined gt -900.0),nindx)
 
    x = hii[indx].zstrong_12oh_oiiinii_pettini
    xerr = hii[indx].zstrong_12oh_oiiinii_pettini_err

    y = hii[indx].zstrong_12oh_kd02_combined
    yerr = hii[indx].zstrong_12oh_kd02_combined_err

    ytitle = '12 + log (O/H) KD02'

    stats = im_stats(x-y,/verbose)
    xstr = strtrim(string(stats.median,format='(F12.2)'),2)+'\pm'+strtrim(string(stats.sigma_rej,format='(F12.2)'),2)
    
    zindicators_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, legendtype=0, charsize=charsize_6, $
      xstyle=3, ystyle=3, /right, /top, position=pos[*,2], /noerase
    djs_oplot, !x.crange, !y.crange, thick=postthick
 
    legend, textoidl(xstr), /right, /bottom, box=0, charsize=charsize_6, charthick=postthick
    legend, '(c)', /left, /top, box=0, charsize=charsize_6, charthick=postthick
    
; the title

;   xyouts, pos[0,0]*0.25, 0.5*(pos[3,0]-pos[1,2])+pos[1,2], textoidl(ytitle), $
;     orientation=90, align=0.5, charsize=charsize_6, charthick=postthick, /normal
;   xyouts, 0.5*(pos[2,1]-pos[0,0])+pos[0,0], pos[1,2]*0.15, textoidl(xtitle), align=0.5, $
;     charsize=charsize_6, charthick=postthick, /normal

    im_openclose, postscript=postscript, /close    

return
end    
