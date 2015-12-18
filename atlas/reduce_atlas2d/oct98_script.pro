pro oct98_script, doplot=doplot, zptshift=zptshift, web=web
; jm02jan3uofa

    dateroot = '98oct'
    dates = dateroot+['15','16','17','18','19']

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
    calibfile = 'caliblist_'+dates+'.txt'
    tellfile = 'stdlist_'+dates+'.txt'
    tellfits = 'telluric_'+dates+'.fits'
    
; ---------------------------------------------------------------------------
; pre-processing
; ---------------------------------------------------------------------------
    
;   iheader_check, root='a.'
;   ispec_makelists, ['a.11','a.12','a.13','a.14','a.15'], $
;     ['98oct15','98oct16','98oct17','98oct18','98oct19'], /all, /overwrite, /gzip
    
;   skyflatfile = 'skyflat_98oct.fits'
;   skyflats = file_search('a.'+['110[3-4]','120[3-4]','130[3-4]','140[3-4]','150[4-5]']+'.fits')
;   arm_flatcombine, skyflats, skyflatfile, [1208,10,1217,132], /display, $
;     psfile='qaplot_skyflat_98oct.ps'
    
;   imake_keywordfile, ['SCANLEN','PA','APERTURE'], comments=['scan length [arcsec]',$
;     'slit position angle [degrees]','slit aperture [arcsec]'], values=[0.0,90.0,2.5], $
;     root='a.01', headfile=headfile

; correct headers

;   h = headfits('a.1133.fits')
;   sxaddpar, h, 'OBJECT', 'BD+284211 2.5'
;   modfits, 'a.1133.fits', 0, h
;   
;   h = headfits('a.1107.fits')
;   sxaddpar, h, 'OBJECT', 'Wolf1346 4.5'
;   modfits, 'a.1107.fits', 0, h
;   
;   h = headfits('a.1224.fits')
;   sxaddpar, h, 'OBJECT', 'NGC 695 drift'
;   modfits, 'a.1224.fits', 0, h
    
; update the headers
    
;   headfile = 'header_keywords_98oct.dat'
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

; NGC0157 drift - no shift needed
    icrcombine, 'ra.'+['1124','1125']+'.fits', /find_pixshift, $
      user_pixshift=[0.0,0.0], /wfits;, /debug
    ibatch, paramfile[0], crlist='ra.1124_2.fits', /crclean
    
; NGC7742 drift - no shift needed
    icrcombine, 'ra.'+['1321','1322']+'.fits', /find_pixshift, $
      user_pixshift=[0.0,0.0], /wfits;, /debug
    ibatch, paramfile[0], crlist='ra.1321_2.fits', /crclean

; IRAS17208-0014 drift
    icrcombine, 'ra.'+['1512','1513']+'.fits', /find_pixshift, $
      user_pixshift=[0.0,1.0], /wfits;, /debug
    ibatch, paramfile[0], crlist='ra.1512_2.fits', /crclean

; repair LA_COSMIC pixels

    repair_la_cosmic_crpix, ['cra.1135_2.fits','cra.1325_2.fits'], $
      ['repair.1135_2.dat','repair.1325_2.dat'], /repair

; ---------------------------------------------------------------------------    
; select sky apertures and sky subtract
; ---------------------------------------------------------------------------    

;   iskyselect, skyfile=skyfile, skyapfile=skyapfile, norder_sky=1, /skyremember

    iall, paramfile, skyapfile=skyapfile

; ---------------------------------------------------------------------------    
; make the sensitivity curves
; ---------------------------------------------------------------------------    

; generate sensitivity curves for each night (2.5")

    senslist = ['sens_2.5_98oct15.fits','sens_2.5_98oct16.fits','sens_2.5_98oct17.fits',$
      'sens_2.5_98oct18.fits','sens_2.5_98oct19.fits']
    senstitle = '1998 '+['October 15','October 16','October 17','October 18','October 19']+' (2.5" Slit)'
    iall, paramfile, stdfile=stdfile2_5, senstitle=senstitle, sensinfo=sensinfo, grey=2, $
      slit_width=2.5, doplot=doplot, sensname=senslist

;   readcol, stdfile2_5[i], stdlist, format='A', /silent, comment='#'
;   info = isensfunc(i1dnames(stdlist,aperture=20.5),grey=2,sensname=senslist[i],$
;     senstitle=senstitle[i],doplot=doplot,/write,slit_width=2.5,zptshift=0.0)

; generate sensitivity curves for each night (4.5")

    senstitle = '1998 '+['October 15','October 16','October 17','October 18','October 19']
    iall, paramfile, stdfile=stdfile4_5, senstitle=senstitle, sensinfo=sensinfo, grey=2, $
      slit_width=4.5, doplot=doplot

;   senslist = ['sens_4.5_98oct15.fits','sens_4.5_98oct16.fits','sens_4.5_98oct17.fits',$
;     'sens_4.5_98oct18.fits','sens_4.5_98oct19.fits']
;   readcol, stdfile[i], stdlist, format='A', /silent, comment='#'
;   info = isensfunc(i1dnames(stdlist,aperture=20.5),grey=2,sensname=senslist[i],$
;     senstitle=senstitle[i],doplot=doplot,/write,slit_width=4.5,zptshift=0.0)

; generate a sensitivity curve for the drift-scanned stars

    ibatch, paramfile[3], stdlist=['scra.1347.fits','scra.1412.fits'], $
      senstitle='1998 October (Drift)', /makesens, sensname='sens_drift_98oct.fits', $
      grey=2, slit_width=20.0, doplot=doplot

; make the combined sensitivity functions
    
    sensname = 'sens_2.5_98oct.fits' & senstitle='1998 October (2.5" Slit)'
    readcol, stdallfile2_5, stdlist, format='A', /silent, comment='#'
    ibatch, paramfile[0], stdlist=stdlist, sensname=sensname, /makesens, zptshift=0.0, $
      grey=2, sensinfo=sensinfo, senstitle=senstitle, slit_width=2.5, doplot=doplot, $
      wmapname='wmaplist_98oct.txt'

;   info = isensfunc(i1dnames(stdlist,aperture=20.5),grey=2,sensname=sensname,$
;     senstitle=senstitle,doplot=doplot,/write,slit_width=2.5,zptshift=0.0)

    sensname = 'sens_4.5_98oct.fits' & senstitle='1998 October (4.5" Slit)'
    readcol, stdallfile4_5, stdlist, format='A', /silent, comment='#'
    ibatch, paramfile[0], stdlist=stdlist, sensname=sensname, /makesens, zptshift=zptshift, $
      grey=2, sensinfo=sensinfo, senstitle=senstitle, slit_width=4.5, doplot=doplot, $
      wmapname='wmaplist_98oct.txt'

;   info = isensfunc(i1dnames(stdlist,aperture=20.5),grey=2,sensname=sensname,$
;     senstitle=senstitle,doplot=doplot,/write,slit_width=4.5,zptshift=zptshift)

; ---------------------------------------------------------------------------    
; generate the telluric spectra
; ---------------------------------------------------------------------------    

    for i = 0L, n_elements(paramfile)-1L do begin

       readcol, tellfile[i], stdlist, format='A', /silent, comment='#'
       hotlist = i1dnames(stdlist,aperture=20.5)

       iconstruct_telluric, hotlist, tellmethod=2L, contmethod=1L, method=3L, $
         npoly=2, psname='qaplot_'+repstr(tellfits[i],'.fits','.ps'), tellfits=tellfits[i], $
         doplot=doplot, debug=0, /write

    endfor
    
; ---------------------------------------------------------------------------    
; flux-calibrate
; ---------------------------------------------------------------------------    

    sensname = 'sens_4.5_98oct.fits'
    iall, paramfile, fluxcalibfile=calibfile, calibfile=calibfile, $
      sensname=sensname;, tellfits=tellfits

; ---------------------------------------------------------------------------    
; miscellaneous tasks
; ---------------------------------------------------------------------------    

; compare the sensitivity curves (ISENSCOMPARE needs to be generalized
; to handle variable-sized sensitivity functions)

    label = [ ['15 October','Clear'],['16 October','Clear'],['17 October','Clear'],$
      ['18 October','Clear'],['19 October','Late clouds'],['Drift',''],['2.5" Slit',''] ]
    senslist = ['sens_4.5_98oct15.fits','sens_4.5_98oct16.fits','sens_4.5_98oct17.fits',$
      'sens_4.5_98oct18.fits','sens_4.5_98oct19.fits','sens_drift_98oct.fits','sens_2.5_98oct.fits']
    meansens = 'sens_4.5_98oct.fits'
    isenscompare, senslist, meansens=meansens, title='1998 October', $
      label=label, psname='qaplot_sens_compare_98oct.ps', /postscript

    loglist = ['qalog_sens_4.5_98oct15.log','qalog_sens_4.5_98oct16.log',$
      'qalog_sens_4.5_98oct17.log','qalog_sens_4.5_98oct18.log','qalog_sens_4.5_98oct19.log',$
      'qalog_sens_drift_98oct.log','qalog_sens_2.5_98oct.log','qalog_sens_4.5_98oct.log']
    label = [ ['15 October','Clear'],['16 October','Clear'],['17 October','Clear'],$
      ['18 October','Clear'],['19 October','Late clouds'],['Drift',''],['2.5" Slit',''],['4.5" Slit',''] ]
    iqaphoto, loglist, label=label, title='1998 October', psname='qaplot_photo_98oct.ps', /postscript
    
; ---------------------------------------------------------------------------    
; generate the web page
; ---------------------------------------------------------------------------    

    if keyword_set(web) then begin
       readcol, 'weblist_98oct.txt', weblist, format='A', comment='#', /silent
       ispec_webpage, '98oct', weblist=weblist, html_path=atlas_path(/dataweb), $
         html_only=html_only, /cleanpng
    endif

; ---------------------------------------------------------------------------    
; call the initial reductions using IBATCH
; ---------------------------------------------------------------------------    

;   ibatch, paramfile[0], /makeflat, checkillum=checkillum, doplot=doplot
;   
;   readcol, procfile[0], proclist, format='A', /silent, comment='#'
;   ibatch, paramfile[0], proclist=proclist, /ccdproc, checkoverscan=checkoverscan
;   
;;  readcol, tracefile[0], tracelist, format='A', /silent, comment='#'
;;  ibatch, paramfile[0], tracelist=tracelist, /distortion, doplot=doplot
;
;   ibatch, paramfile[0], /arcfit, doplot=doplot
;   
;   ibatch, paramfile[0], skyapfile=skyapfile[0], /skysub, doplot=doplot
;   
;   ibatch, paramfile[0], crsplitfile=crsplitfile[0], /crsplits
;   
;   readcol, crfile[0], crlist, format='A', /silent, comment='#'
;   ibatch, paramfile[0], crlist=crlist, /crclean

; ---------------------------------------------------------------------------    
; old code
; ---------------------------------------------------------------------------    

; fix specific images for cosmic rays

;   ibatch, paramfile[1], crlist='ra.1219_2.fits', /crclean, objlim=1.5 ; MRK347
;   ibatch, paramfile[2], crlist='ra.1319_2.fits', /crclean, objlim=1.8 ; N7800
;   ibatch, paramfile[2], crlist='ra.1325_2.fits', /crclean, objlim=2.3 ; IC1623
;   ibatch, paramfile[2], crlist='ra.1333_2.fits', /crclean, objlim=1.0

return
end
