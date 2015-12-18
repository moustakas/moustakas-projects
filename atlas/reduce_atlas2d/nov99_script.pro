pro nov99_script, doplot=doplot, zptshift=zptshift, web=web
; jm02jan3uofa

    dateroot = '99nov'
    dates = dateroot+['02','03','04','29']

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
    
;  iheader_check, root='a.'
;  ispec_makelists, ['a.24','a.25','a.26','a.27'], ['99nov02','99nov03','99nov04','99nov29'], $
;    /all, /overwrite, /gzip

;   skyflats = file_search('a.'+['240[4-6]','250[4-5]','260[3-4]']+'.fits')
;   skyflatfile = 'skyflat_99nov.fits'
;   arm_flatcombine, skyflats, skyflatfile, [1208,5,1217,130], /display, psfile='qaplot_skyflat_99nov.ps'
    
;   imake_keywordfile, ['SCANLEN','POSANGLE','APERTURE'], comments=['scan length [arcsec]',$
;     'slit position angle [degrees]','slit aperture [arcsec]'], values=[0.0,90.0,2.5], $
;     root='a.01', headfile=headfile

; correct headers

;   h = headfits('a.2620.fits')
;   sxaddpar, h, 'OBJECT', 'Feige110 4.5'
;   modfits, 'a.2620.fits', 0, h
    
; repair a.2606.fits for bad rows

;   image = float(readfits('a.2606.fits',header,/silent))
;   mask = byte(image*0.0)
;   mask[151:346,93:96] = 1B
;   cleanimage = djs_maskinterp(image,mask,iaxis=1L)
;   sxaddhist, 'Pixels [151:346,93:96] repaired '+im_today()+'.', header
;   writefits, 'a.2606.fits', cleanimage, header

; update the headers

;   headfile = 'header_keywords_99nov.dat'
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

; ESO602-025 drift; CR on top of OIII5007(?)
    icrcombine, 'ra.'+['2517','2518']+'.fits', /find_pixshift, $
      user_pixshift=[0.0,-1.0], /wfits, nsig=[5,3,3,3,2];, /debug
    ibatch, paramfile[1], crlist='ra.2517_2.fits', /crclean

; UGC12150 drift
    icrcombine, 'ra.'+['2611','2612']+'.fits', /find_pixshift, $
      user_pixshift=[0.0,0.0], /wfits;, /debug
    ibatch, paramfile[2], crlist='ra.2611_2.fits', /crclean

; UGC02238 drift
    icrcombine, 'ra.'+['2617','2618']+'.fits', /find_pixshift, $
      user_pixshift=[0.0,1.0], /wfits;, /debug
    ibatch, paramfile[2], crlist='ra.2617_2.fits', /crclean

; clean up residual cosmic rays

    ibatch, paramfile[0], crlist='ra.2415.fits', /crclean, sigclip=10.0, objlim=1.5, iaxis=1
    ibatch, paramfile[1], crlist='ra.2513.fits', /crclean, sigclip=10.0, objlim=1.0, iaxis=1
    ibatch, paramfile[1], crlist='ra.2514.fits', /crclean, sigclip=10.0, objlim=1.0, iaxis=1
    ibatch, paramfile[1], crlist='ra.2517_2.fits', /crclean, sigclip=5.0, objlim=1.0, iaxis=1

; repair LA_COSMIC pixels

    repair_la_cosmic_crpix, ['cra.2415.fits','cra.2517_2.fits','cra.2611_2.fits',$
      'cra.2622_2.fits','cra.2627_2.fits','cra.2629_2.fits'], ['repair.2415.dat',$
      'repair.2517_2.dat','repair.2611_2.dat','repair.2622_2.dat','repair.2627_2.dat',$
      'repair.2629_2.dat'], /repair

; ---------------------------------------------------------------------------    
; select sky apertures and sky subtract
; ---------------------------------------------------------------------------    

;   iskyselect, skyfile=skyfile, skyapfile=skyapfile, norder_sky=1, /skyremember

    iall, paramfile, skyapfile=skyapfile

; ---------------------------------------------------------------------------    
; make the sensitivity curves
; ---------------------------------------------------------------------------    

; generate sensitivity curves for each night (2.5")

    senslist = ['sens_2.5_99nov02.fits','sens_2.5_99nov03.fits','sens_2.5_99nov04.fits','sens_2.5_99nov29.fits']
    senstitle = '1999 '+['November 02','November 03','November 04','November 29']+' (2.5" Slit)'
    iall, paramfile, stdfile=stdfile2_5, senstitle=senstitle, sensinfo=sensinfo, $
      grey=2, slit_width=2.5, doplot=doplot, sensname=senslist

;   readcol, stdfile2_5[i], stdlist, format='A', /silent, comment='#'
;   info = isensfunc(i1dnames(stdlist,aperture=21),grey=2,sensname=senslist[i],$
;     senstitle=senstitle[i],doplot=doplot,/write,slit_width=2.5,zptshift=0.0)

; generate sensitivity curves for each night except November 29
; (4.5"); November 29 must be treated separately because of ZPTSHIFT

    senstitle = '1999 '+['November 02','November 03','November 04','']
    iall, paramfile, stdfile=stdfile4_5, senstitle=senstitle, sensinfo=sensinfo, $
      grey=2, slit_width=4.5, doplot=doplot, rednight=[0,1,2]

;   senslist = ['sens_4.5_99nov02.fits','sens_4.5_99nov03.fits','sens_4.5_99nov04.fits','sens_4.5_99nov29.fits']
;   readcol, stdfile[i], stdlist, format='A', /silent, comment='#'
;   info = isensfunc(i1dnames(stdlist,aperture=21),grey=2,sensname=senslist[i],$
;     senstitle=senstitle[i],doplot=doplot,/write,slit_width=4.5,zptshift=0.0)

; generate the sensitivity curve2 for November 29, including ZPTSHIFT 

    senstitle = '1999 '+['','','','November 29']
    iall, paramfile, stdfile=stdfile4_5, senstitle=senstitle, sensinfo=sensinfo, $
      grey=2, slit_width=4.5, zptshift=zptshift, doplot=doplot, rednight=3
    
; generate a sensitivity curve for the drift-scanned star

    ibatch, paramfile[2], stdlist='scra.2633.fits', sensname='sens_drift_99nov.fits', $
      /makesens, grey=2, senstitle='1999 November (Drift)', slit_width=20.0, $
      doplot=doplot
    
; make the combined sensitivity function

    sensname = 'sens_2.5_99nov.fits' & senstitle = '1999 November (2.5" Slit)'
    readcol, stdallfile2_5, stdlist, format='A', /silent, comment='#'
    ibatch, paramfile[0], stdlist=stdlist, sensname=sensname, /makesens, zptshift=0.0, $
      grey=2, sensinfo=sensinfo, senstitle=senstitle, slit_width=2.5, doplot=doplot, $
      wmapname='wmaplist_99nov.txt'

;   info = isensfunc(i1dnames(stdlist,aperture=21),grey=2,sensname=sensname,$
;     senstitle=senstitle,doplot=doplot,/write,slit_width=2.5,zptshift=0.0)

    sensname = 'sens_4.5_99nov.fits' & senstitle='1999 November (4.5" Slit)'
    readcol, stdallfile4_5, stdlist, format='A', /silent, comment='#'
    ibatch, paramfile[0], stdlist=stdlist, sensname=sensname, /makesens, zptshift=zptshift, $
      grey=2, sensinfo=sensinfo, senstitle=senstitle, slit_width=4.5, doplot=doplot, $
      wmapname='wmaplist_99nov.txt'

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

    sensname = 'sens_4.5_99nov.fits'
    iall, paramfile, fluxcalibfile=calibfile, calibfile=calibfile, $
      sensname=sensname, rednight=[0,1,2];, tellfits=tellfits

    sensname = 'sens_4.5_99nov29.fits'
    iall, paramfile, fluxcalibfile=calibfile, calibfile=calibfile, $
      sensname=sensname, rednight=3;, tellfits=tellfits

; ---------------------------------------------------------------------------    
; miscellaneous tasks
; ---------------------------------------------------------------------------    

; compare the sensitivity curves (for Nov 2-4 + drift only)

    label = [ ['02 November','Thin clouds'],['03 November','Clear'],['04 November','Clear'],$
      ['Drift',''],['2.5" Slit',''] ]
    senslist = ['sens_4.5_99nov02.fits','sens_4.5_99nov03.fits',$
      'sens_4.5_99nov04.fits','sens_drift_99nov.fits','sens_2.5_99nov.fits']
    meansens = 'sens_4.5_99nov.fits'
    isenscompare, senslist, meansens=meansens, title='1999 November', $ 
      label=label, psname='qaplot_sens_compare_99nov.ps', /postscript

    loglist = ['qalog_sens_4.5_99nov02.log','qalog_sens_4.5_99nov03.log',$
      'qalog_sens_4.5_99nov04.log','qalog_sens_drift_99nov.log',$
      'qalog_sens_2.5_99nov.log','qalog_sens_4.5_99nov.log']
    label = [ ['02 November','Thin clouds'],['03 November','Clear'],['04 November','Clear'],$
      ['Drift',''],['2.5" Slit',''],['4.5" Slit',''] ]
    iqaphoto, loglist, label=label, title='1999 November', psname='qaplot_photo_99nov.ps', /postscript

    loglist = ['qalog_sens_2.5_99nov29.log','qalog_sens_4.5_99nov29.log']
    label = [ ['29 November','2.5" Slit','Cirrus'],['29 November','4.5" Slit','Cirrus'] ]
    iqaphoto, loglist, label=label, title='1999 November 29', psname='qaplot_photo_99nov29.ps', /postscript
    
; ---------------------------------------------------------------------------    
; generate the web page
; ---------------------------------------------------------------------------    

    if keyword_set(web) then begin
       readcol, 'weblist_99nov.txt', weblist, format='A', comment='#', /silent
       ispec_webpage, '99nov', weblist=weblist, html_path=atlas_path(/dataweb), $
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

;   ibatch, paramfile[0], /crclean, crlist='ra.2415.fits', sigfrac=1.5 ; N7713
;   ibatch, paramfile[0], /crclean, crlist='ra.2416.fits', sigfrac=1.5 ; N7640
;   ibatch, paramfile[0], /crclean, crlist='ra.2422.fits', sigfrac=1.5 ; N784
;   ibatch, paramfile[0], /crclean, crlist='ra.2424.fits', objlim=1.2 ; N959
;   ibatch, paramfile[1], /crclean, crlist='ra.2517_2.fits', sigfrac=1.0 ; ESO602-025
;   ibatch, paramfile[1], /crclean, crlist='ra.2521.fits', objlim=1.6, sigfrac=1.0 ; UGC12588
;   ibatch, paramfile[1], /crclean, crlist='ra.2523_2.fits', objlim=1.6 ; IC5298=ZW475-056
;   icrcombine, ['ra.2531.fits','ra.2532.fits'], nsig=[3,3,2,2,1.5], /wfits ; U04881
;   ibatch, paramfile[1], /crclean, crlist='ra.2531_2.fits'

;   ibatch, paramfile[2], crlist='ra.2611_2.fits', /crclean, sigfrac=1.5 ; U12150
;   ibatch, paramfile[2], crlist='ra.2614.fits', /crclean, objlim=1.5 ; U01561
;   ibatch, paramfile[2], crlist='ra.2615_2.fits', /crclean, sigfrac=1.0 ; N232
;   ibatch, paramfile[2], crlist='ra.2617_2.fits', /crclean, objlim=1.5 ; U02238
;   ibatch, paramfile[2], crlist='ra.2622_2.fits', /crclean, objlim=1.0, sigclip=4.0 ; N1560
;   ibatch, paramfile[2], crlist='ra.2627_2.fits', /crclean, objlim=1.2 ; N2090
;   ibatch, paramfile[2], crlist='ra.2629_2.fits', /crclean, sigfrac=1.0 ; MCG+08-18-012
;   ibatch, paramfile[3], crlist='ra.2714_2.fits', /crclean, objlim=1.2 ; U00903
;   ibatch, paramfile[3], crlist='ra.2716.fits', /crclean, iaxis=1 ; N615
    
return
end
