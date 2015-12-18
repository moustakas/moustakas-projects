pro may03_script, doplot=doplot, zptshift=zptshift, web=web

    dateroot = '03may'
    dates = dateroot+['26','27']

    paramfile = 'ibatch_'+dates+'.txt'
    procfile = 'objlist_'+dates+'.txt'
    skyfile = 'skylist_'+dates+'.txt'
    skyapfile = 'skyaplist_'+dates+'.txt'
    tracefile = 'tracelist_'+dates+'.txt'
    crsplitfile = 'crsplits_'+dates+'.txt'
    crfile = 'crlist_'+dates+'.txt'
    stdfile4_5 = 'stdlist_4.5_'+dates+'.txt'
    stdfile2_5 = 'stdlist_2.5_'+dates+'.txt'
    stdallfile4_5 = 'stdlist_4.5_'+dateroot+'.txt'
    stdallfile2_5 = 'stdlist_2.5_'+dateroot+'.txt'
    stddriftfile = 'stdlist_drift_'+dateroot+'.txt'
    calibfile = 'caliblist_'+dates+'.txt'
    tellfile = 'stdlist_'+dates+'.txt'
    tellfits = 'telluric_'+dates+'.fits'
    
; ---------------------------------------------------------------------------
; pre-processing
; ---------------------------------------------------------------------------
    
;   iheader_check, root='a.'
;   ispec_makelists, ['a.50','a.51'], ['03may26','03may27'], /all, /overwrite, /gzip
    
;   skyflats = file_search('a.'+['500[3-5]','510[1-3]']+'.fits')
;   skyflatfile = 'skyflat_03may.fits'
;   arm_flatcombine, skyflats, skyflatfile, [1209,15,1218,273], /display, psfile='qaplot_skyflat_03may.ps'
    
; May 26
;   
;   biasfile = 'bias_03may26.fits'
;   biases = file_search('a.5090[0-2]?.fits')
;   arm_zerocombine, biases, biasfile, /display, psfile='qaplot_bias_03may26.ps'
;
;   domefile = 'domeflat_03may26.fits'
;   domeflats = file_search('a.5090[5-6]?.fits')
;   arm_flatcombine, domeflats, domefile, [1205,0,1217,269], /display, psfile='qaplot_domeflat_03may26.ps'
;
; May 27
;   
;   biasfile = 'bias_03may27.fits'
;   biases = file_search('a.5190[0-2]?.fits')
;   arm_zerocombine, biases, biasfile, /display, psfile='qaplot_bias_03may27.ps'
;
;   domefile = 'domeflat_03may27.fits'
;   domeflats = file_search('a.5190[5-6]?.fits')
;   arm_flatcombine, domeflats, domefile, [1209,15,1218,273], /display, psfile='qaplot_domeflat_03may27.ps'
    
;   imake_keywordfile, ['SCANLEN','PA','APERTURE'], comments=['scan length [arcsec]',$
;     'slit position angle [degrees]','slit aperture [arcsec]'], values=[0.0,90.0,2.5], $
;     root='a.01', headfile=headfile

;   headfile = 'header_keywords_03may.dat'
;   iheader_keywords, headfile, /update, silent=silent

; ---------------------------------------------------------------------------
; carry out all the initial reductions
; ---------------------------------------------------------------------------
    
    iall, paramfile, procfile=procfile, doplot=doplot, tracefile=tracefile

; ---------------------------------------------------------------------------    
; combine cr-splits and cosmic-ray reject
; ---------------------------------------------------------------------------    

;   ibatch, paramfile[0], crsplitfile=crsplitfile[0], /find_pixshift, /debug, /crsplits

    iall, paramfile, crsplitfile=crsplitfile, crfile=crfile, /find_pixshift

; MRK0331 drift 
    icrcombine, 'ra.'+['5039','5040']+'.fits', /find_pixshift, $
      user_pixshift=[0.0,-4.0], /wfits;, /debug
    ibatch, paramfile[0], crlist='ra.5039_2.fits', /crclean

; ---------------------------------------------------------------------------    
; select sky apertures and sky subtract
; ---------------------------------------------------------------------------    

;   iskyselect, skyfile=skyfile, skyapfile=skyapfile, norder_sky=1, /skyremember

    iall, paramfile, skyapfile=skyapfile

; ---------------------------------------------------------------------------    
; make the sensitivity curves
; ---------------------------------------------------------------------------    

; generate sensitivity curves for each night (2.5")

    senslist = ['sens_2.5_03may26.fits','sens_2.5_03may27.fits']
    senstitle = '2003 '+['May 26','May 27']+' (2.5" Slit)'
    iall, paramfile, stdfile=stdfile2_5, senstitle=senstitle, sensinfo=sensinfo, $
      slit_width=2.5, grey=2, doplot=doplot, sensname=senslist

;   readcol, stdfile2_5[i], stdlist, format='A', /silent, comment='#'
;   info = isensfunc(i1dnames(stdlist,aperture=21.58),grey=2,sensname=senslist[i],$
;     senstitle=senstitle[i],doplot=doplot,/write,slit_width=4.5,zptshift=0.0)

; generate sensitivity curves for each night (4.5")

    senstitle = '2003 '+['May 26','May 27']
    iall, paramfile, stdfile=stdfile4_5, senstitle=senstitle, sensinfo=sensinfo, $
      slit_width=4.5, grey=2, doplot=doplot

;   senslist = ['sens_4.5_03may26.fits','sens_4.5_03may27.fits']
;   readcol, stdfile[i], stdlist, format='A', /silent, comment='#'
;   info = isensfunc(i1dnames(stdlist,aperture=21.58),grey=2,sensname=senslist[i],$
;     senstitle=senstitle[i],doplot=doplot,/write,slit_width=4.5,zptshift=0.0)

; generate the drift-scanned sensitivity function    

    sensname = 'sens_drift_03may.fits' & senstitle='2003 May (Drift)'
    readcol, stddriftfile, stdlist, format='A', /silent, comment='#'
    ibatch, paramfile[0], stdlist=stdlist, senstitle=senstitle, $
      sensname=sensname, /makesens, grey=2, slit_width=20.0, doplot=doplot, $
      wmapname='wmaplist_03may.txt'

;   info = isensfunc(i1dnames(stdlist,aperture=21.58),grey=2,sensname=sensname,$
;     senstitle=senstitle,doplot=doplot,/write,slit_width=20,zptshift=0.0)

; make the combined sensitivity functions
    
    sensname = 'sens_2.5_03may.fits' & senstitle='2003 May (2.5" Slit)'
    readcol, stdallfile2_5, stdlist, format='A', /silent, comment='#'
    ibatch, paramfile[0], stdlist=stdlist, sensname=sensname, /makesens, zptshift=0.0, $
      grey=2, sensinfo=sensinfo, senstitle=senstitle, slit_width=2.5, doplot=doplot, $
      wmapname='wmaplist_03may.txt'
    
;   info = isensfunc(i1dnames(stdlist,aperture=21.58),grey=2,sensname=sensname,$
;     senstitle=senstitle,doplot=doplot,/write,slit_width=2.5,zptshift=0.0)

    sensname = 'sens_4.5_03may.fits' & senstitle='2003 May (4.5" Slit)'
    readcol, stdallfile4_5, stdlist, format='A', /silent, comment='#'
    ibatch, paramfile[0], stdlist=stdlist, sensname=sensname, /makesens, zptshift=zptshift, $
      grey=2, sensinfo=sensinfo, senstitle=senstitle, slit_width=4.5, doplot=doplot, $
      wmapname='wmaplist_03may.txt'

;   info = isensfunc(i1dnames(stdlist,aperture=21.58),grey=2,sensname=sensname,$
;     senstitle=senstitle,doplot=doplot,/write,slit_width=4.5,zptshift=zptshift)

; ---------------------------------------------------------------------------    
; generate the telluric spectra
; ---------------------------------------------------------------------------    

    for i = 0L, n_elements(paramfile)-1L do begin

       readcol, tellfile[i], stdlist, format='A', /silent, comment='#'
       hotlist = i1dnames(stdlist,aperture=21.58)

       iconstruct_telluric, hotlist, tellmethod=2L, contmethod=1L, method=3L, $
         npoly=2, psname='qaplot_'+repstr(tellfits[i],'.fits','.ps'), tellfits=tellfits[i], $
         doplot=doplot, debug=0, /write

    endfor
    
; ---------------------------------------------------------------------------    
; flux-calibrate
; ---------------------------------------------------------------------------    
    
    sensname = 'sens_4.5_03may.fits'
    iall, paramfile, fluxcalibfile=calibfile, calibfile=calibfile, $
      sensname=sensname;, tellfits=tellfits

; ---------------------------------------------------------------------------    
; miscellaneous tasks
; ---------------------------------------------------------------------------    

; compare the sensitivity curves

    label = [['26 May','Clear'],['27 May','Early Clouds'],['Drift',''],['2.5" Slit',''] ]
    senslist = ['sens_4.5_03may26.fits','sens_4.5_03may27.fits',$
      'sens_drift_03may.fits','sens_2.5_03may.fits']
    meansens = 'sens_4.5_03may.fits'
    isenscompare, senslist, meansens=meansens, label=label, title='2003 May', $
      psname='qaplot_sens_compare_03may.ps', /postscript

    loglist = ['qalog_sens_4.5_03may26.log','qalog_sens_4.5_03may27.log','qalog_sens_drift_03may.log',$
      'qalog_sens_2.5_03may.log','qalog_sens_4.5_03may.log']
    label = [['26 May','Clear'],['27 May','Early Clouds'],['Drift',''],['2.5" Slit',''],['4.5" slit','']]
    iqaphoto, loglist, label=label, title='2003 May', psname='qaplot_photo_03may.ps', /postscript
    
; ---------------------------------------------------------------------------    
; generate the web page
; ---------------------------------------------------------------------------    

    if keyword_set(web) then begin
       readcol, 'weblist_03may.txt', weblist, format='A', comment='#', /silent
       ispec_webpage, '03may', weblist=weblist, html_path=atlas_path(/dataweb), $
         html_only=html_only, /cleanpng
    endif

; ---------------------------------------------------------------------------    
;; measure and write out the spectral resolution
; ---------------------------------------------------------------------------    
;
;    readcol, 'arclist_03may26.txt', arclist, format='A', comment='#'
;    ibatch, paramfile[0], caliblist=arclist, sensname='', extfile='', tracename='', /calibrate
;    lampres = im_specres_lamp('w'+arclist,psname='lampres_03may26.ps',/postscript,/write)
;
;    fwhm2sig = 2.0*sqrt(2.0*alog(2.0))
;
;    narc = n_elements(arclist)
;    medwave = djs_median(lampres.linewave,2)
;    medres = fwhm2sig*djs_median(lampres.linewidth,2)
;    
;    wave = findgen(3500)+3500
;    coeff = poly_fit(lampres.linewave,fwhm2sig*lampres.linewidth,2)
;;   coeff = poly_fit(medwave,medres,2)
;    lampfit = poly(wave,float(coeff))
;
;    dfpsplot, 'lampres_median.ps', /square
;    
;    plotsym, 8, 2.0
;    djs_plot, [0], [0], xrange=[3700,6900], yrange=[8,13], /nodata, $
;      xthick=5.0, ythick=5.0, charsize=2.0, charthick=5.0, xsty=3, ysty=3, $
;      xtitle='Wavelength ['+angstrom()+']', ytitle='FWHM Resolution ['+angstrom()+']'
;    for i = 0L, narc-1L do djs_oplot, lampres[i].linewave, fwhm2sig*lampres[i].linewidth, ps=8
;;   djs_oplot, medwave, medres, ps=8
;    djs_oplot, wave, lampfit, line=0, thick=3.0
;
;    dfpsclose
;
;    out = {specwave: wave, specres: lampfit}
;
;    outname = atlas_path(/specfit)+'atlas1d_specres.fits'
;    mwrfits, out, outname, /create
        
return 
end   
