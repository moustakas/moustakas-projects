pro deep2_oiilf, vmax_zbin1, vmax_zbin2, vmax_zbin3, vmax_zbin4, $
  vmax_willmer, postscript=postscript
; jm07sep30nyu - written
; jm08may02nyu - major update    

; read the data

    deepparent = read_deep2_oiilf_sample(/parent)
    deepkcorr = read_deep2_oiilf_sample(/kcorr)
    deepidust = read_deep2_oiilf_sample(/ispec)

    sdssparent = read_sdss_oiilf_sample(/parent)
    sdsskcorr = read_sdss_oiilf_sample(/kcorr)
    sdssidust = read_sdss_oiilf_sample(/ispec)
    
; compute the spectroscopic solid angle    
    
;   splog, 'Computing DEEP2/DR3 area.'
;   im_deep2_area, dr3_area
    dr3_area = 0.00085281834 ; 3.0*!dtor^2.0

; constants

    red, h100=0.7, omega0=0.3, omega_lambda=0.7
    h100 = redh100()
    lsun = 3.826D33
    
    timelabel1 = [1.0,3.0,5.0,7.0,9] ; [Gyr]

; initialize some path names and plotting variables

    oiilfpath = deep2_path(/projects)+'oiilf/'
    paperpath = deep2_path(/papers)+'oiilf/'
    pspath = paperpath+'FIG_OIILF/'

    postthick1 = 2.0
    postthick2 = 2.0
    postthick3 = 2.0
    postthick4 = 2.0
    postthick5 = 2.0
    postthick6 = 1.0
    postthick7 = 2.0
    
    loadct, 0, /silent
    textcolor1 = 'white'
    axis_color = 'white'
    scolor = 'forest green'
    dcolor = 'dodger blue'
    
    if (not keyword_set(postscript)) and (not keyword_set(pdf)) then $
      im_window, 0, xratio=0.5, /square

    if keyword_set(pdf) then begin

       postscript = 0L
       encapsulated = 0L
       pspath = paperpath+'/keynote/'

       loadct, 3, /silent
       textcolor1 = 'white'
       axis_color = 'black'

       postthick1 = 6.0 
       postthick2 = 6.0 
       postthick3 = 8.0 
       postthick4 = 4.0 
       postthick5 = 8.0 
       postthick6 = 1.0 
       postthick7 = 10.0 

       agespsize = 0.1 & agessym = 108 & agescolor = 'black'
       sdsspsize = 0.2 & sdsssym = 108 & sdsscolor = 'black' 

       scolor = 'forest green'
       dcolor = 'dodger blue'

    endif
    
    if keyword_set(postscript) then begin
    
       loadct, 0, /silent
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

    charsize_0 = 1.0
    charsize_1 = 1.1
    charsize_2 = 1.2
    charsize_3 = 1.3
    charsize_4 = 1.4
    charsize_5 = 1.5
    charsize_6 = 1.6
    charsize_7 = 1.7
    charsize_8 = 1.8
    charsize_9 = 1.9
    singlecharsize_0 = 2.0
    singlecharsize_1 = 2.1
    singlecharsize_2 = 2.2
    singlecharsize_3 = 2.3
    singlecharsize_4 = 2.4
    singlecharsize_5 = 2.5
    charsize_30 = 3.0

    phirange1 = alog10([1D-6,1D-2]) ; alog10([1D-5,3D-2])
    phirange2 = alog10([1D-6,5D-2]) ; alog10([1D-5,3D-2])
    mbrange1 = [-19.0,-24.5]
    mbrange2 = [-18.5,-24.5]
    loiirange1 = [39.3,44.7]
    loiirange2 = [40.1,43.4]
    ewoiirange1 = [-0.5,3.5]
    loiibinsize = 0.2 ; 0.05
    rrange1 = [20.0,25.0]
    mbbinsize = 0.25
    
; some sample definitions    

    loiilocut_zbin1 = 40.5  & loiihicut_zbin1 = 43.0
    loiilocut_zbin2 = 40.75 & loiihicut_zbin2 = 43.0
    loiilocut_zbin3 = 41.0  & loiihicut_zbin3 = 43.0
    loiilocut_zbin4 = 41.25 & loiihicut_zbin4 = 43.0
    
    mblocut_zbin1 = -19.5 & mbhicut_zbin1 = -24.0
    mblocut_zbin2 = -20.0 & mbhicut_zbin2 = -24.0
    mblocut_zbin3 = -20.5 & mbhicut_zbin3 = -24.0
    mblocut_zbin4 = -21.0 & mbhicut_zbin4 = -24.0
    
    lozcut_zbin1 = 0.7 & hizcut_zbin1 = 0.9
    lozcut_zbin2 = 0.9 & hizcut_zbin2 = 1.1
    lozcut_zbin3 = 1.1 & hizcut_zbin3 = 1.3
    lozcut_zbin4 = 1.3 & hizcut_zbin4 = 1.5
    
;;; read the VMAX values in each bin
;;
;;    splog, 'Reading the VMAX files for each redshift bin.'
;;    
;;    if (n_elements(vmax_zbin1) eq 0L) or (n_elements(vmax_zbin2) eq 0L) or $
;;      (n_elements(vmax_zbin3) eq 0L) or (n_elements(vmax_zbin4) eq 0L) then begin
;;
;;       vmax_zbin1 = mrdfits(analysis_path+'deep2_vmax_0.7_0.9.fits',1,/silent)
;;       spherematch, deep2dust.ra, deep2dust.dec, vmax_zbin1.ra, vmax_zbin1.dec, $
;;         0.1/3600.0, junk, vmax_indx, maxmatch=1
;;
;;       vmax_zbin1 = (vmax_zbin1)[vmax_indx]
;;       vmax_zbin2 = (mrdfits(analysis_path+'deep2_vmax_0.9_1.1.fits',1,/silent))[vmax_indx]
;;       vmax_zbin3 = (mrdfits(analysis_path+'deep2_vmax_1.1_1.3.fits',1,/silent))[vmax_indx]
;;       vmax_zbin4 = (mrdfits(analysis_path+'deep2_vmax_1.3_1.5.fits',1,/silent))[vmax_indx]
;;
;;    endif
;;
;;    if (n_elements(vmax_willmer) eq 0L) then begin
;;
;;       vmax_willmer = mrdfits(analysis_path+'deep2_vmax_0.8_1.0.fits',1,/silent)
;;       spherematch, deep2parent.ra, deep2parent.dec, vmax_willmer.ra, vmax_willmer.dec, $
;;         0.1/3600.0, junk, vmax_bigindx, maxmatch=1
;;       vmax_willmer = vmax_willmer[vmax_bigindx]
;;
;;    endif

;   vmax_zbin1 = vmax_willmer
;   vmax_zbin2 = vmax_willmer
;   vmax_zbin3 = vmax_willmer
;   vmax_zbin4 = vmax_willmer
    
; ---------------------------------------------------------------------------
; SDSS [O II] LF - we are plotting (1/V)*dN/d(log L)
; ---------------------------------------------------------------------------

    psname = 'sdss_lf_loii'
    xpage = 8.5 & ypage = 8.5
    im_openclose, pspath+psname, postscript=(keyword_set(postscript) or $
      keyword_set(pdf)), xsize=xpage, ysize=ypage, encapsulated=encapsulated, $
      silent=keyword_set(pdf)

    pagemaker, nx=1, ny=1, yspace=0.0, xspace=0.0, width=7.0, height=7.0, $
      xmargin=[1.2,0.3], ymargin=[0.4,1.1], xpage=8.5, ypage=8.5, position=pos, /normal

    xrange = [5.5,10.0]+alog10(lsun)
    yrange = [-9.0,-1.0]

    xtitle = textoidl('log L([O II]) (erg s^{-1})')
    ytitle = textoidl('log \Phi'+'(log L[O II]) (Mpc^{-3})')
    charsize = 1.4

    sdssindx = where(sdsskcorr.vmax_noevol gt 0.0,nsdssindx)
    x = sdssidust[sdssindx].oii_lum+alog10(lsun)
    weight = 1.0/sdsskcorr[sdssindx].vmax_noevol

    ybin = im_hist1d(x,weight,binsize=loiibinsize,$
      obin=xbin,binedge=0,h_err=ybin_err,/ylog)
    ybin = ybin - alog10(loiibinsize)

    djs_plot, [0], [0], /nodata, xthick=postthick1, ythick=postthick1, $
      xtitle=xtitle, ytitle=ytitle, charsize=charsize_9, charthick=postthick2, $
      xsty=1, ysty=1, xrange=xrange, position=pos[*,0], yrange=yrange
    oploterror, xbin, ybin, ybin_err, ps=10, line=0, thick=postthick1, errthick=postthick1

stop    
    
;   djs_oplot, xbin1, ybin1, ps=10, line=0, thick=postthick1
    legend, '(a) 0.7 < z < 0.9', /left, /top, box=0, charsize=1.5, charthick=postthick2

    im_openclose, postscript=(keyword_set(postscript) or keyword_set(pdf)), /close    
    if keyword_set(pdf) then begin
       splog, 'Writing '+pspath+psname+'.pdf'
       spawn, 'ps2pdf13 '+pspath+psname+'.ps '+pspath+psname+'.pdf', /sh
       rmfile, pspath+psname+'.ps' ; delete the junk ps file
    endif

stop    
    
; ---------------------------------------------------------------------------
; EW([O II]) vs z 
; ---------------------------------------------------------------------------

    psname = 'deep_z_vs_ewoii'
    xpage = 8.5 & ypage = 8.5
    im_openclose, pspath+psname, postscript=(keyword_set(postscript) or $
      keyword_set(pdf)), xsize=xpage, ysize=ypage, encapsulated=encapsulated, $
      silent=keyword_set(pdf)

    pagemaker, nx=1, ny=1, yspace=0.0, xspace=0.0, width=7.0, height=7.0, $
      xmargin=[1.2,0.3], ymargin=[0.4,1.1], xpage=8.5, ypage=8.5, position=pos, /normal

    xrange = [0.7,1.5]
    yrange = ewoiirange1

    xtitle = 'Redshift'
    ytitle = textoidl('log EW([O II]) (\AA)')
    charsize = 2.0

    indx = where((deepkcorr.z gt 0.7) and (deepkcorr.z lt 1.5),nzbin1)
    x = deepkcorr[indx].z
    y = alog10(deepidust[indx].oii_3727_ew[0])

    im_symbols, 108, psize=0.4, /fill
    djs_plot, [0], [0], /nodata, xthick=postthick1, ythick=postthick1, $
      xtitle=xtitle, ytitle=ytitle, charsize=charsize, charthick=postthick2, $
      xsty=3, ysty=1, xrange=xrange, position=pos[*,0], yrange=yrange
    djs_oplot, x, y, ps=8

    im_openclose, postscript=(keyword_set(postscript) or keyword_set(pdf)), /close    
    if keyword_set(pdf) then begin
       splog, 'Writing '+pspath+psname+'.pdf'
       spawn, 'ps2pdf13 '+pspath+psname+'.ps '+pspath+psname+'.pdf', /sh
       rmfile, pspath+psname+'.ps' ; delete the junk ps file
    endif

stop    
    
; ---------------------------------------------------------------------------
; EW([O II]) vs R
; ---------------------------------------------------------------------------

    psname = 'deep_R_vs_ewoii'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.5, encapsulated=encapsulated

    pagemaker, nx=1, ny=1, yspace=0.0, xspace=0.0, width=7.0, height=7.0, $
      xmargin=[1.2,0.3], ymargin=[0.4,1.1], xpage=8.5, ypage=8.5, position=pos, /normal

    xrange = rrange1
    yrange = ewoiirange1

    xtitle = 'R (AB mag)'
    ytitle = textoidl('log EW([O II]) (\AA)')
    charsize = 2.0

    indx = where((deepkcorr.z gt 0.7) and (deepkcorr.z lt 1.5),nzbin1)
    x = deepkcorr[indx].magr
    y = alog10(deepidust[indx].oii_3727_ew[0])

    djs_plot, [0], [0], /nodata, xthick=postthick1, ythick=postthick1, $
      xtitle=xtitle, ytitle=ytitle, charsize=charsize, charthick=postthick2, $
      xsty=3, ysty=1, xrange=xrange, position=pos[*,0], yrange=yrange
    djs_oplot, x, y, ps=8

    im_openclose, postscript=postscript, /close    

; ---------------------------------------------------------------------------
; bivariate LF in four redshift bins
; ---------------------------------------------------------------------------

    psname = 'deep_loii_vs_mb'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.5, encapsulated=encapsulated

    pagemaker, nx=2, ny=2, yspace=0, xspace=0, width=3.5*[1,1], height=3.5*[1,1], $
      xmargin=[1.3,0.2], ymargin=[0.4,1.1], xpage=8.5, ypage=8.5, position=pos, /normal

    xrange = loiirange1
    yrange = mbrange2
    xtitle = textoidl('log L([O II]) (erg s^{-1})')
    ytitle = textoidl('M_{B} (mag)')
    charsize = 1.7
    
; -------------------------    
; 0.7<z<0.9
; -------------------------    

    zbin1 = where((deepkcorr.z gt 0.7) and (deepkcorr.z lt 0.9) and (vmax_zbin1.vmax gt 0.0),nzbin1)
    x1 = deepidust[zbin1].oii_lum
    y1 = deepkcorr[zbin1].ubvri_absmag[1]
    weight1 = 1.0/((dr3_area/3.0)*vmax_zbin1[zbin1].vmax)
    
    hogg_scatterplot, x1, y1, weight=weight1, xthick=postthick1, ythick=postthick1, $
      xtitle='', xtickname=replicate(' ',10), $
      ytitle=ytitle, charsize=charsize, charthick=postthick2, xrange=xrange, $
      yrange=yrange, xsty=1, ysty=1, position=pos[*,0], /outliers, outpsym=8, outsymsize=1.0
    legend, '(a) 0.7 < z < 0.9', /left, /top, box=0, charsize=1.5, charthick=postthick2
; ZBIN1
;   djs_oplot, [loiilocut_zbin1,loiihicut_zbin1], mblocut_zbin1*[1,1], line=0, thick=postthick3
;   djs_oplot, [loiilocut_zbin1,loiihicut_zbin1], mbhicut_zbin1*[1,1], line=0, thick=postthick3
;   djs_oplot, loiilocut_zbin1*[1,1], [mblocut_zbin1,mbhicut_zbin1], line=0, thick=postthick3
;   djs_oplot, loiihicut_zbin1*[1,1], [mblocut_zbin1,mbhicut_zbin1], line=0, thick=postthick3

; -------------------------    
; 0.9<z<1.1
; -------------------------    

    zbin2 = where((deepkcorr.z gt 0.9) and (deepkcorr.z lt 1.1) and (vmax_zbin2.vmax gt 0.0),nzbin2)
    x2 = deepidust[zbin2].oii_lum
    y2 = deepkcorr[zbin2].ubvri_absmag[1]
    weight2 = 1.0/((dr3_area/3.0)*vmax_zbin2[zbin2].vmax)
;   weight2 = 1.0/deepkcorr[zbin2].weight
    
    hogg_scatterplot, x2, y2, weight=weight2, /overplot, /noerase, $
      xthick=postthick1, ythick=postthick1, xtitle='', xtickname=replicate(' ',10), $
      ytitle='', ytickname=replicate(' ',10), charsize=charsize, charthick=postthick2, xrange=xrange, $
      yrange=yrange, xsty=1, ysty=1, position=pos[*,1], /outliers, outpsym=8, outsymsize=1.0
    legend, '(b) 0.9 < z < 1.1', /left, /top, box=0, charsize=1.5, charthick=postthick2

; -------------------------    
; 1.1<z<1.3
; -------------------------    

    zbin3 = where((deepkcorr.z gt 1.1) and (deepkcorr.z lt 1.3) and (vmax_zbin3.vmax gt 0.0),nzbin3)
    x3 = deepidust[zbin3].oii_lum
    y3 = deepkcorr[zbin3].ubvri_absmag[1]
    weight3 = 1.0/((dr3_area/3.0)*vmax_zbin3[zbin3].vmax)
;   weight3 = 1.0/deepkcorr[zbin3].weight
    
    hogg_scatterplot, x3, y3, weight=weight3, /overplot, /noerase, xthick=postthick1, ythick=postthick1, xtitle=xtitle, $
      ytitle=ytitle, charsize=charsize, charthick=postthick2, xrange=xrange, $
      yrange=yrange, xsty=1, ysty=1, position=pos[*,2], /outliers, outpsym=8, outsymsize=1.0
    legend, '(c) 1.1 < z < 1.3', /left, /top, box=0, charsize=1.5, charthick=postthick2

; -------------------------    
; 1.3<z<1.5
; -------------------------    

    zbin4 = where((deepkcorr.z gt 1.3) and (deepkcorr.z lt 1.5) and (vmax_zbin4.vmax gt 0.0),nzbin4)
    x4 = deepidust[zbin4].oii_lum
    y4 = deepkcorr[zbin4].ubvri_absmag[1]
    weight4 = 1.0/((dr3_area/3.0)*vmax_zbin4[zbin4].vmax)
;   weight4 = 1.0/deepkcorr[zbin4].weight
    
    hogg_scatterplot, x4, y4, weight=weight4, /overplot, /noerase, xthick=postthick1, ythick=postthick1, xtitle=xtitle, $
      ytitle='', ytickname=replicate(' ',10), charsize=charsize, charthick=postthick2, xrange=xrange, $
      yrange=yrange, xsty=1, ysty=1, position=pos[*,3], /outliers, outpsym=8, outsymsize=1.0
    legend, '(d) 1.3 < z < 1.5', /left, /top, box=0, charsize=1.5, charthick=postthick2

    im_openclose, postscript=postscript, /close    

; ---------------------------------------------------------------------------
; L([O II]) vs z 
; ---------------------------------------------------------------------------

    psname = 'deep_z_vs_loii'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.5, encapsulated=encapsulated

    pagemaker, nx=1, ny=1, yspace=0.0, xspace=0.0, width=7.0, height=7.0, $
      xmargin=[1.2,0.3], ymargin=[0.4,1.1], xpage=8.5, ypage=8.5, position=pos, /normal

    xrange = [0.7,1.5]
    yrange = loiirange2

    xtitle = 'Redshift'
    ytitle = textoidl('log L([O II]) (erg s^{-1})')
    charsize = 2.0

    indx = where((deepkcorr.z gt 0.7) and (deepkcorr.z lt 1.5),nzbin1)
    x = deepkcorr[indx].z
    y = deepidust[indx].oii_lum

    djs_plot, [0], [0], /nodata, xthick=postthick1, ythick=postthick1, $
      xtitle=xtitle, ytitle=ytitle, charsize=charsize, charthick=postthick2, $
      xsty=3, ysty=1, xrange=xrange, position=pos[*,0], yrange=yrange
    djs_oplot, x, y, ps=8
    factor = 1.005
; ZBIN1
    djs_oplot, [lozcut_zbin1,hizcut_zbin1], loiilocut_zbin1*[1,1], line=0, thick=postthick3
    djs_oplot, [lozcut_zbin1,hizcut_zbin1], loiihicut_zbin1*[1,1], line=0, thick=postthick3
    djs_oplot, lozcut_zbin1*[1,1], [loiilocut_zbin1,loiihicut_zbin1], line=0, thick=postthick3
    djs_oplot, hizcut_zbin1*[1,1], [loiilocut_zbin1,loiihicut_zbin1], line=0, thick=postthick3
; ZBIN2
    djs_oplot, [lozcut_zbin2,hizcut_zbin2], loiilocut_zbin2*[1,1], line=2, thick=postthick3
    djs_oplot, [lozcut_zbin2,hizcut_zbin2], loiihicut_zbin2*[1,1], line=2, thick=postthick3
    djs_oplot, lozcut_zbin2*[1,1]*factor, [loiilocut_zbin2,loiihicut_zbin2], line=2, thick=postthick3
    djs_oplot, hizcut_zbin2*[1,1], [loiilocut_zbin2,loiihicut_zbin2], line=2, thick=postthick3
; ZBIN3
    djs_oplot, [lozcut_zbin3,hizcut_zbin3], loiilocut_zbin3*[1,1], line=3, thick=postthick3
    djs_oplot, [lozcut_zbin3,hizcut_zbin3], loiihicut_zbin3*[1,1], line=3, thick=postthick3
    djs_oplot, lozcut_zbin3*[1,1]*factor, [loiilocut_zbin3,loiihicut_zbin3], line=3, thick=postthick3
    djs_oplot, hizcut_zbin3*[1,1], [loiilocut_zbin3,loiihicut_zbin3], line=3, thick=postthick3
; ZBIN4
    djs_oplot, [lozcut_zbin4,hizcut_zbin4], loiilocut_zbin4*[1,1], line=1, thick=postthick3
    djs_oplot, [lozcut_zbin4,hizcut_zbin4], loiihicut_zbin4*[1,1], line=1, thick=postthick3
    djs_oplot, lozcut_zbin4*[1,1]*factor, [loiilocut_zbin4,loiihicut_zbin4], line=1, thick=postthick3
    djs_oplot, hizcut_zbin4*[1,1], [loiilocut_zbin4,loiihicut_zbin4], line=1, thick=postthick3

    im_openclose, postscript=postscript, /close    

; ---------------------------------------------------------------------------
; [O II] LF in four redshift bins; we are plotting (1/V)*dN/d(log L)
; ---------------------------------------------------------------------------

    psname = 'deep_lf_loii'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.5, encapsulated=encapsulated

    pagemaker, nx=2, ny=2, yspace=0, xspace=0, width=3.5*[1,1], height=3.5*[1,1], $
      xmargin=[1.1,0.4], ymargin=[0.4,1.1], xpage=8.5, ypage=8.5, position=pos, /normal

    xrange = loiirange2
    yrange = phirange2

    xtitle = textoidl('log L([O II]) (erg s^{-1})')
    ytitle = textoidl('log \Phi'+'(log L[O II]) (Mpc^{-3})')
    charsize = 1.4

; -------------------------    
; 0.7<z<0.9
; -------------------------    

    zbin1 = where((deepkcorr.z gt 0.7) and (deepkcorr.z lt 0.9) and (vmax_zbin1.vmax gt 0.0) and $
      (deepkcorr.ubvri_absmag[1] lt mblocut_zbin1) and (deepkcorr.ubvri_absmag[1] gt mbhicut_zbin1) and $
      (deepidust.oii_lum gt loiilocut_zbin1) and (deepidust.oii_lum lt loiihicut_zbin1),nzbin1)

    vol1 = (dr3_area/3.0)*((lf_comvol(0.9)-lf_comvol(0.7))[0])/h100^3.0 ; h=0.7
    x1 = deepidust[zbin1].oii_lum
    weight1 = 1.0/((dr3_area/3.0)*vmax_zbin1[zbin1].vmax)

    ybin1 = im_hist1d(x1,weight1,binsize=loiibinsize,obin=xbin1,binedge=0,$
      h_err=ybin1_err,/ylog)
    ybin1 = ybin1 - alog10(loiibinsize)

    djs_plot, [0], [0], /nodata, xthick=postthick1, ythick=postthick1, $
      xtitle='', ytitle=ytitle, charsize=charsize, charthick=postthick2, $
      xsty=1, ysty=1, xrange=xrange, position=pos[*,0], yrange=yrange, $
      xtickname=replicate(' ',10)
    oploterror, xbin1, ybin1, ybin1_err, ps=10, line=0, thick=postthick1, errthick=postthick1
;   djs_oplot, xbin1, ybin1, ps=10, line=0, thick=postthick1
    legend, '(a) 0.7 < z < 0.9', /left, /top, box=0, charsize=1.5, charthick=postthick2

; -------------------------    
; 0.9<z<1.1
; -------------------------    

    zbin2 = where((deepkcorr.z gt 0.9) and (deepkcorr.z lt 1.1) and (vmax_zbin2.vmax gt 0.0) and $
      (deepkcorr.ubvri_absmag[1] lt mblocut_zbin2) and (deepkcorr.ubvri_absmag[1] gt mbhicut_zbin2) and $
      (deepidust.oii_lum gt loiilocut_zbin2) and (deepidust.oii_lum lt loiihicut_zbin2),nzbin2)

    vol2 = (dr3_area/3.0)*((lf_comvol(1.1)-lf_comvol(0.9))[0])/h100^3.0
    x2 = deepidust[zbin2].oii_lum
    weight2 = 1.0/((dr3_area/3.0)*vmax_zbin2[zbin2].vmax)

    ybin2 = im_hist1d(x2,weight2,binsize=loiibinsize,obin=xbin2,binedge=0,$
      h_err=ybin2_err,/ylog)
    ybin2 = ybin2 - alog10(loiibinsize)

    djs_plot, [0], [0], /nodata, /noerase, xthick=postthick1, ythick=postthick1, $
      xtitle='', ytitle='', charsize=charsize, charthick=postthick2, $
      xsty=1, ysty=1, xrange=xrange, position=pos[*,1], yrange=yrange, $
      xtickname=replicate(' ',10), ytickname=replicate(' ',10)
    djs_oplot, xbin1, ybin1, ps=10, line=1, thick=postthick1
    oploterror, xbin2, ybin2, ybin2_err, ps=10, line=0, thick=postthick1, errthick=postthick1
    legend, '(b) 0.9 < z < 1.1', /left, /top, box=0, charsize=1.5, charthick=postthick2

; -------------------------    
; 1.1<z<1.3
; -------------------------    

    zbin3 = where((deepkcorr.z gt 1.1) and (deepkcorr.z lt 1.3) and (vmax_zbin3.vmax gt 0.0) and $
      (deepkcorr.ubvri_absmag[1] lt mblocut_zbin3) and (deepkcorr.ubvri_absmag[1] gt mbhicut_zbin3) and $
      (deepidust.oii_lum gt loiilocut_zbin3) and (deepidust.oii_lum lt loiihicut_zbin3),nzbin3)

    vol3 = (dr3_area/3.0)*((lf_comvol(1.3)-lf_comvol(1.1))[0])/h100^3.0
    x3 = deepidust[zbin3].oii_lum
    weight3 = 1.0/((dr3_area/3.0)*vmax_zbin3[zbin3].vmax)

    ybin3 = im_hist1d(x3,weight3,binsize=loiibinsize,obin=xbin3,binedge=0,$
      h_err=ybin3_err,/ylog)
    ybin3 = ybin3 - alog10(loiibinsize)

    djs_plot, [0], [0], /nodata, /noerase, xthick=postthick1, ythick=postthick1, $
      xtitle=xtitle, ytitle=ytitle, charsize=charsize, charthick=postthick2, $
      xsty=1, ysty=1, xrange=xrange, yrange=yrange, position=pos[*,2];, xtickinterval=1.0
    djs_oplot, xbin1, ybin1, ps=10, line=1, thick=postthick1
    oploterror, xbin3, ybin3, ybin3_err, ps=10, line=0, thick=postthick1, errthick=postthick1
;   djs_oplot, xbin3, ybin3, ps=10, line=0, thick=postthick1
    legend, '(c) 1.1 < z < 1.3', /left, /top, box=0, charsize=1.5, charthick=postthick2

; -------------------------    
; 1.3<z<1.5
; -------------------------    

    zbin4 = where((deepkcorr.z gt 1.3) and (deepkcorr.z lt 1.5) and (vmax_zbin4.vmax gt 0.0) and $
      (deepkcorr.ubvri_absmag[1] lt mblocut_zbin4) and (deepkcorr.ubvri_absmag[1] gt mbhicut_zbin4) and $
      (deepidust.oii_lum gt loiilocut_zbin4) and (deepidust.oii_lum lt loiihicut_zbin4),nzbin4)

    vol4 = (dr3_area/3.0)*((lf_comvol(1.5)-lf_comvol(1.3))[0])/h100^3.0
    x4 = deepidust[zbin4].oii_lum
    weight4 = 1.0/((dr3_area/3.0)*vmax_zbin4[zbin4].vmax)

    ybin4 = im_hist1d(x4,weight4,binsize=loiibinsize,obin=xbin4,binedge=0,$
      h_err=ybin4_err,/ylog)
    ybin4 = ybin4 - alog10(loiibinsize)

    djs_plot, [0], [0], /nodata, /noerase, xthick=postthick1, ythick=postthick1, $
      xtitle=xtitle, xrange=xrange, ytitle='', charsize=charsize, charthick=postthick2, $
      xsty=1, ysty=1, yrange=yrange, position=pos[*,3], ytickname=replicate(' ',10);, xtickinterval=1.0
    djs_oplot, xbin1, ybin1, ps=10, line=1, thick=postthick1
    oploterror, xbin4, ybin4, ybin4_err, ps=10, line=0, thick=postthick1, errthick=postthick1
;   djs_oplot, xbin4, ybin4, ps=10, line=0, thick=postthick1
    legend, '(d) 1.3 < z < 1.5', /left, /top, box=0, charsize=1.5, charthick=postthick2

    im_openclose, postscript=postscript, /close    

; ---------------------------------------------------------------------------
; M_B LF in four redshift bins; we are plotting (1/V)*dN/d(log L)
; ---------------------------------------------------------------------------

    psname = 'deep_lf_mb'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.5, encapsulated=encapsulated

    pagemaker, nx=2, ny=2, yspace=0, xspace=0, width=3.5*[1,1], height=3.5*[1,1], $
      xmargin=[1.1,0.4], ymargin=[0.4,1.1], xpage=8.5, ypage=8.5, position=pos, /normal

    xrange = mbrange1
    yrange = phirange2

    xtitle = textoidl('M_{B} (mag)')
    ytitle = textoidl('log \Phi'+'(M_{B}) (Mpc^{-3} mag^{-1})')
    charsize = 1.7

; -------------------------    
; 0.7<z<0.9
; -------------------------    

    zbin1 = where((deepkcorr.z gt 0.7) and (deepkcorr.z lt 0.9) and (vmax_zbin1.vmax gt 0.0) and $
      (deepkcorr.ubvri_absmag[1] lt mblocut_zbin1) and (deepkcorr.ubvri_absmag[1] gt mbhicut_zbin1),nzbin1)

    vol1 = (dr3_area/3.0)*((lf_comvol(0.9)-lf_comvol(0.7))[0])/h100^3.0 ; h=1-->h=0.7
    x1 = deepkcorr[zbin1].ubvri_absmag[1]
    weight1 = 1.0/((dr3_area/3.0)*vmax_zbin1[zbin1].vmax)

    ybin1 = im_hist1d(x1,weight1,binsize=mbbinsize,obin=xbin1,binedge=0,$
      histmax=mblocut_zbin1,histmin=mbhicut_zbin1,h_err=ybin1_err,/ylog)
    ybin1 = ybin1 - alog10(mbbinsize)

    djs_plot, [0], [0], /nodata, xthick=postthick1, ythick=postthick1, $
      xtitle='', ytitle=ytitle, charsize=charsize, charthick=postthick2, $
      xsty=1, ysty=1, xrange=xrange, position=pos[*,0], yrange=yrange, $
      xtickname=replicate(' ',10);;, /ylog
    oploterror, xbin1, ybin1, ybin1_err, ps=10, line=0, thick=postthick1, errthick=postthick1
;   djs_oplot, xbin1, ybin1, ps=10, line=0, thick=postthick1
    legend, '(a) 0.7 < z < 0.9', /left, /top, box=0, charsize=1.5, charthick=postthick2

; -------------------------    
; 0.9<z<1.1
; -------------------------    

    zbin2 = where((deepkcorr.z gt 0.9) and (deepkcorr.z lt 1.1) and (vmax_zbin2.vmax gt 0.0) and $
      (deepkcorr.ubvri_absmag[1] lt mblocut_zbin2) and (deepkcorr.ubvri_absmag[1] gt mbhicut_zbin2),nzbin2)

    vol2 = (dr3_area/3.0)*((lf_comvol(1.1)-lf_comvol(0.9))[0])/h100^3.0
    x2 = deepkcorr[zbin2].ubvri_absmag[1]
    weight2 = 1.0/((dr3_area/3.0)*vmax_zbin2[zbin2].vmax)

    ybin2 = im_hist1d(x2,weight2,binsize=mbbinsize,obin=xbin2,binedge=0,$
      histmax=mblocut_zbin2,histmin=mbhicut_zbin2,h_err=ybin2_err,/ylog)
    ybin2 = ybin2 - alog10(mbbinsize)

    djs_plot, [0], [0], /nodata, /noerase, xthick=postthick1, ythick=postthick1, $
      xtitle='', ytitle='', charsize=charsize, charthick=postthick2, $
      xsty=1, ysty=1, xrange=xrange, position=pos[*,1], yrange=yrange, $
      xtickname=replicate(' ',10), ytickname=replicate(' ',10);, /ylog
    djs_oplot, xbin1, ybin1, ps=10, line=1, thick=postthick1
    oploterror, xbin2, ybin2, ybin2_err, ps=10, line=0, thick=postthick1, errthick=postthick1
;   djs_oplot, xbin2, ybin2, ps=10, line=0, thick=postthick1
    legend, '(b) 0.9 < z < 1.1', /left, /top, box=0, charsize=1.5, charthick=postthick2

; -------------------------    
; 1.1<z<1.3
; -------------------------    

    zbin3 = where((deepkcorr.z gt 1.1) and (deepkcorr.z lt 1.3) and (vmax_zbin3.vmax gt 0.0) and $
      (deepkcorr.ubvri_absmag[1] lt mblocut_zbin3) and (deepkcorr.ubvri_absmag[1] gt mbhicut_zbin3),nzbin3)

    vol3 = (dr3_area/3.0)*((lf_comvol(1.3)-lf_comvol(1.1))[0])/h100^3.0
    x3 = deepkcorr[zbin3].ubvri_absmag[1]
    weight3 = 1.0/((dr3_area/3.0)*vmax_zbin3[zbin3].vmax)

    ybin3 = im_hist1d(x3,weight3,binsize=mbbinsize,obin=xbin3,binedge=0,$
      histmax=mblocut_zbin3,histmin=mbhicut_zbin3,h_err=ybin3_err,/ylog)
    ybin3 = ybin3 - alog10(mbbinsize)

    djs_plot, [0], [0], /nodata, /noerase, xthick=postthick1, ythick=postthick1, $
      xtitle=xtitle, ytitle=ytitle, charsize=charsize, charthick=postthick2, $
      xsty=1, ysty=1, xrange=xrange, yrange=yrange, position=pos[*,2];;, /ylog
    djs_oplot, xbin1, ybin1, ps=10, line=1, thick=postthick1
    oploterror, xbin3, ybin3, ybin3_err, ps=10, line=0, thick=postthick1, errthick=postthick1
;   djs_oplot, xbin3, ybin3, ps=10, line=0, thick=postthick1
    legend, '(c) 1.1 < z < 1.3', /left, /top, box=0, charsize=1.5, charthick=postthick2

; -------------------------    
; 1.3<z<1.5
; -------------------------    

    zbin4 = where((deepkcorr.z gt 1.3) and (deepkcorr.z lt 1.5) and (vmax_zbin4.vmax gt 0.0) and $
      (deepkcorr.ubvri_absmag[1] lt mblocut_zbin4) and (deepkcorr.ubvri_absmag[1] gt mbhicut_zbin4),nzbin4)

    vol4 = (dr3_area/3.0)*((lf_comvol(1.5)-lf_comvol(1.3))[0])/h100^3.0
    x4 = deepkcorr[zbin4].ubvri_absmag[1]
    weight4 = 1.0/((dr3_area/3.0)*vmax_zbin4[zbin4].vmax)

    ybin4 = im_hist1d(x4,weight4,binsize=mbbinsize,obin=xbin4,binedge=0,$
      histmax=mblocut_zbin4,histmin=mbhicut_zbin4,h_err=ybin4_err,/ylog)
    ybin4 = ybin4 - alog10(mbbinsize)

    djs_plot, [0], [0], /nodata, /noerase, xthick=postthick1, ythick=postthick1, $
      xtitle=xtitle, xrange=xrange, ytitle='', charsize=charsize, charthick=postthick2, $
      xsty=1, ysty=1, yrange=yrange, position=pos[*,3], ytickname=replicate(' ',10);;, /ylog
    djs_oplot, xbin1, ybin1, ps=10, line=1, thick=postthick1
    oploterror, xbin4, ybin4, ybin4_err, ps=10, line=0, thick=postthick1, errthick=postthick1
;   djs_oplot, xbin4, ybin4, ps=10, line=0, thick=postthick1
    legend, '(d) 1.3 < z < 1.5', /left, /top, box=0, charsize=1.5, charthick=postthick2

    im_openclose, postscript=postscript, /close    

; ---------------------------------------------------------------------------
; compare my M_B LF at 0.8<z<1.0 with Willmer et al. 2006
; ---------------------------------------------------------------------------

    psname = 'deep_lf_mb_vs_willmer'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.5, encapsulated=encapsulated

    pagemaker, nx=1, ny=1, yspace=0.0, xspace=0.0, width=7.0, height=7.0, $
      xmargin=[1.2,0.3], ymargin=[0.4,1.1], xpage=8.5, ypage=8.5, position=pos, /normal

    xrange = mbrange1
    yrange = phirange1

    xtitle = textoidl('M_{B} (mag)')
    ytitle = textoidl('log \Phi'+'(M_{B}) (Mpc^{-3} mag^{-1})')
    charsize = 2.0

    this_mblocut = -19.5
    this_mbhicut = -23.5
    
    indx = where((deep2parent.z gt 0.8) and (deep2parent.z lt 1.2) and (vmax_willmer.vmax gt 0.0) and $
      (deep2parent.ubvri_absmag[1] lt this_mblocut) and (deep2parent.ubvri_absmag[1] gt this_mbhicut),nindx)

    x1 = deep2parent[indx].ubvri_absmag[1]
    weight1 = 1.0/((dr3_area/3.0)*vmax_willmer[indx].vmax)

    ybin1 = im_hist1d(x1,weight1,binsize=mbbinsize,obin=xbin1,binedge=0,$
      histmax=this_mblocut,histmin=this_mbhicut,h_err=ybin1_err,/ylog)
    ybin1 = ybin1 - alog10(mbbinsize)

    djs_plot, [0], [0], /nodata, xthick=postthick1, ythick=postthick1, $
      xtitle=xtitle, ytitle=ytitle, charsize=charsize, charthick=postthick2, $
      xsty=1, ysty=1, xrange=xrange, position=pos[*,0], yrange=yrange
    oploterror, xbin1, ybin1, ybin1_err, ps=10, line=0, thick=postthick1, errthick=postthick1
    niceprint, xbin1, ybin1, ybin1_err

;   adq_weighted_hist, x1, weight1, -24.5, -19.5, mbbinsize, yy, xx, yyerr
;   djs_oplot, xx, alog10(yy>1D-10), ps=10, color='navy', thick=2
;   adq_weighted_hist, x1, weight1*0.0+1.0, -24.5, -19.5, mbbinsize, yy, xx, yyerr
;   ybin1 = im_hist1d(x1,weight1*0.0+1.0,binsize=mbbinsize,obin=xbin1,binedge=0,histmax=-19.5,histmin=-24.5)
;   djs_plot, xbin1, ybin1, ps=10, xsty=3, ysty=3
;   djs_oplot, xx, yy, ps=10, color='red'
;   plothist, x1, bin=mbbinsize, omax=-19.5, omin=-24.5, color=djs_icolor('dark green'), /overplot, thick=4

; compare with Willmer et al. 2006    
    
    willmer_all_file = getenv('EVOLUTION_DIR')+'/data/paper2/deep2_optimal_all_vmax_0.80_1.00.da'
    readcol, willmer_all_file, mb, logphi, logphi_fit, phi, phi_err, /silent

    djs_oplot, mb, logphi, ps=10, color='dark green', line=5, thick=postthick3

    legend, 'Willmer et al.', line=5, thick=postthick3, /right, /top, box=0, $
      charsize=2.0, charthick=postthick2, color=djs_icolor('dark green');, $
;     position=[xrange[1]+0.2,yrange[1]-0.35]
    legend, '0.8 < z < 1.0', /left, /bottom, box=0, charsize=2.0, $
      charthick=postthick2;, position=[xrange[1]+0.2,yrange[1]-0.15]
    
    im_openclose, postscript=postscript, /close    

; ---------------------------------------------------------------------------
; M_B vs z 
; ---------------------------------------------------------------------------

    psname = 'deep_z_vs_mb'
    im_openclose, pspath+psname, postscript=postscript, xsize=8.5, ysize=8.5, encapsulated=encapsulated

    pagemaker, nx=1, ny=1, yspace=0.0, xspace=0.0, width=7.0, height=7.0, $
      xmargin=[1.2,0.3], ymargin=[0.4,1.1], xpage=8.5, ypage=8.5, position=pos, /normal

    xrange = [0.7,1.5]
    yrange = mbrange2

    xtitle = 'Redshift'
    ytitle = textoidl('M_{B} (mag)')
    charsize = 2.0

    indx = where((deep2parent.z gt 0.7) and (deep2parent.z lt 1.5),nzbin1)
    x = deep2parent[indx].z
    y = deep2parent[indx].ubvri_absmag[1]
;   indx = where((deepkcorr.z gt 0.7) and (deepkcorr.z lt 1.5),nzbin1)
;   x = deepkcorr[indx].z
;   y = deepkcorr[indx].ubvri_absmag[1]

    djs_plot, [0], [0], /nodata, xthick=postthick1, ythick=postthick1, $
      xtitle=xtitle, ytitle=ytitle, charsize=charsize, charthick=postthick2, $
      xsty=3, ysty=1, xrange=xrange, position=pos[*,0], yrange=yrange
    djs_oplot, x, y, ps=8
    factor = 1.005
; ZBIN1
    djs_oplot, [lozcut_zbin1,hizcut_zbin1], mblocut_zbin1*[1,1], line=0, thick=postthick3
    djs_oplot, [lozcut_zbin1,hizcut_zbin1], mbhicut_zbin1*[1,1], line=0, thick=postthick3
    djs_oplot, lozcut_zbin1*[1,1], [mblocut_zbin1,mbhicut_zbin1], line=0, thick=postthick3
    djs_oplot, hizcut_zbin1*[1,1], [mblocut_zbin1,mbhicut_zbin1], line=0, thick=postthick3
; ZBIN2
    djs_oplot, [lozcut_zbin2,hizcut_zbin2], mblocut_zbin2*[1,1], line=2, thick=postthick3
    djs_oplot, [lozcut_zbin2,hizcut_zbin2], mbhicut_zbin2*[1,1], line=2, thick=postthick3
    djs_oplot, lozcut_zbin2*[1,1]*factor, [mblocut_zbin2,mbhicut_zbin2], line=2, thick=postthick3
    djs_oplot, hizcut_zbin2*[1,1], [mblocut_zbin2,mbhicut_zbin2], line=2, thick=postthick3
; ZBIN3
    djs_oplot, [lozcut_zbin3,hizcut_zbin3], mblocut_zbin3*[1,1], line=3, thick=postthick3
    djs_oplot, [lozcut_zbin3,hizcut_zbin3], mbhicut_zbin3*[1,1], line=3, thick=postthick3
    djs_oplot, lozcut_zbin3*[1,1]*factor, [mblocut_zbin3,mbhicut_zbin3], line=3, thick=postthick3
    djs_oplot, hizcut_zbin3*[1,1], [mblocut_zbin3,mbhicut_zbin3], line=3, thick=postthick3
; ZBIN4
    djs_oplot, [lozcut_zbin4,hizcut_zbin4], mblocut_zbin4*[1,1], line=1, thick=postthick3
    djs_oplot, [lozcut_zbin4,hizcut_zbin4], mbhicut_zbin4*[1,1], line=1, thick=postthick3
    djs_oplot, lozcut_zbin4*[1,1]*factor, [mblocut_zbin4,mbhicut_zbin4], line=1, thick=postthick3
    djs_oplot, hizcut_zbin4*[1,1], [mblocut_zbin4,mbhicut_zbin4], line=1, thick=postthick3

;   legend, '(a) 0.7 < z < 0.9', /left, /top, box=0, charsize=1.5, charthick=postthick2

    im_openclose, postscript=postscript, /close    

return
end
