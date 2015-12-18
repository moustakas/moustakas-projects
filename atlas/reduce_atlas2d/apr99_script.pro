pro apr99_script, doplot=doplot, zptshift=zptshift, web=web
; jm02jan3uofa

    dateroot = '99apr'
    dates = dateroot+['19','21','22','23']

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
;   ispec_makelists, ['a.16','a.17','a.18','a.19'], ['99apr19','99apr21','99apr22','99apr23'], $
;     /all, /overwrite, gzip=0
    
;   skyflats = file_search('a.'+['160[4-6]','180[4-5]','190[4-6]']+'.fits')
;   skyflatfile = 'skyflat_99apr.fits'
;   arm_flatcombine, skyflats, skyflatfile, [1208,5,1217,130], /display, psfile='qaplot_skyflat_99apr.ps'
    
;   imake_keywordfile, ['SCANLEN','POSANGLE','APERTURE'], comments=['scan length [arcsec]',$
;     'slit position angle [degrees]','slit aperture [arcsec]'], values=[0.0,90.0,2.5], $
;     root='a.01', headfile=headfile

; fix rows 49 through 51 in night 2 (1999 April 21)

;   rowfix_night2, /wfits

; update the headers    
    
;   headfile = 'header_keywords_99apr.dat'
;   iheader_keywords, headfile, /update, silent=silent

; ---------------------------------------------------------------------------
; carry out all the initial reductions
; ---------------------------------------------------------------------------

    iall, paramfile, procfile=procfile, doplot=doplot, tracefile=tracefile

; ---------------------------------------------------------------------------    
; combine cr-splits and cosmic-ray reject
; ---------------------------------------------------------------------------    

; no pixel shifts are needed since there is only one crsplit file and
; for that object (UGC06436AB) the pixel shift is zero    
    
;   ibatch, paramfile[2], crsplitfile=crsplitfile[2], /find_pixshift, /debug, /crsplits

    iall, paramfile, crsplitfile=crsplitfile, crfile=crfile;, /find_pixshift ; <-- NOTE!

; clean up residual cosmic rays

    ibatch, paramfile[2], crlist='ra.1816.fits', /crclean, objlim=0.5, iaxis=0

; repair LA_COSMIC pixels

    repair_la_cosmic_crpix, 'cra.1629.fits', 'repair.1629.dat', /repair

; ---------------------------------------------------------------------------    
; select sky apertures and sky subtract
; ---------------------------------------------------------------------------    

;   iskyselect, skyfile=skyfile, skyapfile=skyapfile, norder_sky=1, /skyremember
    iall, paramfile, skyapfile=skyapfile

; ---------------------------------------------------------------------------    
; make the sensitivity curves
; ---------------------------------------------------------------------------    

; generate sensitivity curves for all the nights (2.5")

    senstitle = '1999 '+['April 19','April 21','April 22','April 23']+' (2.5" Slit)'
    senslist = ['sens_2.5_99apr19.fits','sens_2.5_99apr21.fits','sens_2.5_99apr22.fits','sens_2.5_99apr23.fits']
    iall, paramfile, stdfile=stdfile2_5, senstitle=senstitle, sensinfo=sensinfo, $
      grey=2, slit_width=2.5, doplot=doplot, sensname=senslist
    
;   readcol, stdfile[i], stdlist, format='A', /silent, comment='#'
;   info = isensfunc(i1dnames(stdlist,aperture=),grey=2,sensname=senslist[i],$
;     senstitle=senstitle[i],doplot=doplot,/write,slit_width=2.5,zptshift=0.0)

; generate sensitivity curves for all the nights (4.5")

    senstitle = '1999 '+['April 19','April 21','April 22','April 23']
    iall, paramfile, stdfile=stdfile4_5, senstitle=senstitle, sensinfo=sensinfo, $
      grey=2, slit_width=4.5, doplot=doplot
    
;   senslist = ['sens_4.5_99apr19.fits','sens_4.5_99apr21.fits','sens_4.5_99apr22.fits','sens_4.5_99apr23.fits']
;   readcol, stdfile[i], stdlist, format='A', /silent, comment='#'
;   info = isensfunc(i1dnames(stdlist,aperture=),grey=2,sensname=senslist[i],$
;     senstitle=senstitle[i],doplot=doplot,/write,slit_width=4.5,zptshift=0.0)

; generate a sensitivity curve for the drift-scanned stars for April

    ibatch, paramfile[1], stdlist=['scra.1717.fits','scra.1911.fits'], $
      sensname='sens_drift_99apr.fits', /makesens, senstitle='1999 April (Drift)', $
      grey=2, sensinfo=sensinfo, slit_width=20.0, doplot=doplot

; make the combined sensitivity function
    
    sensname = 'sens_2.5_99apr.fits' & senstitle='1999 April (2.5" Slit)'
    readcol, stdallfile2_5, stdlist, format='A', /silent, comment='#'
    ibatch, paramfile[0], stdlist=stdlist, sensname=sensname, /makesens, zptshift=0.0, $
      grey=2, sensinfo=sensinfo, senstitle=senstitle, slit_width=2.5, doplot=doplot, $
      wmapname='wmaplist_99apr.txt'

;   info = isensfunc(i1dnames(stdlist,aperture=21),grey=2,sensname=sensname,$
;     senstitle=senstitle,doplot=doplot,/write,slit_width=2.5,zptshift=0.0)

    sensname = 'sens_4.5_99apr.fits' & senstitle='1999 April (4.5" Slit)'
    readcol, stdallfile4_5, stdlist, format='A', /silent, comment='#'
    ibatch, paramfile[0], stdlist=stdlist, sensname=sensname, /makesens, zptshift=zptshift, $
      grey=2, sensinfo=sensinfo, senstitle=senstitle, slit_width=4.5, doplot=doplot, $
      wmapname='wmaplist_99apr.txt'

;   info = isensfunc(i1dnames(stdlist,aperture=21),grey=2,sensname=sensname,$
;     senstitle=senstitle,doplot=doplot,/write,slit_width=4.5,zptshift=zptshift)

; ---------------------------------------------------------------------------    
; generate the telluric spectra
; ---------------------------------------------------------------------------    

    for i = 0L, n_elements(paramfile)-1L do begin

       readcol, tellfile[i], stdlist, format='A', /silent, comment='#'
       hotlist = i1dnames(stdlist,aperture=21)

       iconstruct_telluric, hotlist, tellmethod=2L, contmethod=1L, method=3L, $
         npoly=2, psname='qaplot_'+repstr(tellfits[i],'.fits','.ps'), tellfits=tellfits[i], $
         doplot=doplot, debug=0, /write

    endfor
    
; ---------------------------------------------------------------------------    
; flux-calibrate
; ---------------------------------------------------------------------------    
    
    sensname = 'sens_4.5_99apr.fits'
    iall, paramfile, fluxcalibfile=calibfile, calibfile=calibfile, $
      sensname=sensname;, tellfits=tellfits

; ---------------------------------------------------------------------------    
; miscellaneous tasks
; ---------------------------------------------------------------------------    

; compare the sensitivity curves

    label = [ ['19 April','Clear'],['21 April','Clear, Windy'],['22 April','Cloudy, Windy'],$
      ['23 April','Late clouds'],['Drift',''],['2.5" Slit',''] ]
    senslist = ['sens_4.5_99apr19.fits','sens_4.5_99apr21.fits','sens_4.5_99apr22.fits',$
      'sens_4.5_99apr23.fits','sens_drift_99apr.fits','sens_2.5_99apr.fits']
    meansens = 'sens_4.5_99apr.fits'
    isenscompare, senslist, meansens=meansens, title='1999 April', $ 
      label=label, psname='qaplot_sens_compare_99apr.ps', /postscript

    loglist = ['qalog_sens_4.5_99apr19.log','qalog_sens_4.5_99apr21.log','qalog_sens_4.5_99apr22.log',$
      'qalog_sens_4.5_99apr23.log','qalog_sens_drift_99apr.log','qalog_sens_2.5_99apr.log','qalog_sens_4.5_99apr.log']
    label = [ ['19 April','Clear'],['21 April','Clear, Windy'],['22 April','Cloudy, Windy'],$
      ['23 April','Late clouds'],['Drift',''],['2.5" Slit',''],['4.5" Slit',''] ]
    iqaphoto, loglist, label=label, title='1999 April', psname='qaplot_photo_99apr.ps', /postscript
    
; ---------------------------------------------------------------------------    
; generate the web pages
; ---------------------------------------------------------------------------    

    if keyword_set(web) then begin
       readcol, 'weblist_99apr.txt', weblist, format='A', comment='#', /silent
       ispec_webpage, '99apr', weblist=weblist, html_path=atlas_path(/dataweb), $
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

;   ibatch, paramfile[0], crlist='ra.1617.fits', /crclean, objlim=1.3 ; A1367 pos 7
;   ibatch, paramfile[0], crlist='ra.1624.fits', /crclean, objlim=1.0, iaxis=1 ; N4656
;   ibatch, paramfile[0], crlist='ra.1627.fits', /crclean, objlim=1.5 ; N5430
;   ibatch, paramfile[0], crlist='ra.1629.fits', /crclean, objlim=1.3 ; N5953/4
;   ibatch, paramfile[5], crlist='ra.2111.fits', /crclean, objlim=1.0 ; U06456
;   ibatch, paramfile[5], crlist='ra.2113.fits', /crclean, objlim=1.5 ; N4144
;   ibatch, paramfile[6], crlist='ra.2211.fits', /crclean, objlim=1.2 ; N2500
;   ibatch, paramfile[6], crlist='ra.2215.fits', /crclean, objlim=1.5 ; U07690
;   ibatch, paramfile[6], crlist='ra.2225.fits', /crclean, objlim=1.5 ; IZw107
;   ibatch, paramfile[7], crlist='ra.2311.fits', /crclean, sigfrac=1.7 ; N3432
;   ibatch, paramfile[7], crlist='ra.2313.fits', /crclean, sigfrac=1.7 ; N3104
;   ibatch, paramfile[7], crlist='ra.2315.fits', /crclean, objlim=1.3 ; DDO161
;   ibatch, paramfile[7], crlist='ra.2323.fits', /crclean, objlim=1.3 ; DD0190

return
end
