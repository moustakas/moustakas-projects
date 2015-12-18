pro apr00_script, doplot=doplot, zptshift=zptshift, web=web
; jm02jan3uofa

    dateroot_mar = '00mar'
    dateroot_apr = '00apr'
    dates = [dateroot_mar+['29','30','31'],dateroot_apr+['01','02','03']]

    paramfile = 'ibatch_'+dates+'.txt'
    procfile = 'objlist_'+dates+'.txt'
    skyfile = 'skylist_'+dates+'.txt'
    skyapfile = 'skyaplist_'+dates+'.txt'
    tracefile = 'tracelist_'+dates+'.txt'
    crsplitfile = 'crsplits_'+dates+'.txt'
    crfile = 'crlist_'+dates+'.txt'
    stdfile4_5 = 'stdlist_4.5_'+dates+'.txt'
    stdfile2_5 = 'stdlist_2.5_'+dates+'.txt'
    stdallfile4_5 = 'stdlist_4.5_'+dateroot_apr+'.txt'
    stdallfile2_5 = 'stdlist_2.5_'+dateroot_apr+'.txt'
    calibfile = 'caliblist_'+dates+'.txt'
    tellfile = 'stdlist_'+dates+'.txt'
    tellfits = 'telluric_'+dates+'.fits'

;   dates = '00'+['mar29','mar30','mar31','apr01','apr02','apr03']
;
;   paramfile = 'ibatch_'+dates+'.txt'
;   procfile = 'objlist_'+dates+'.txt'
;   skyfile = 'skylist_'+dates+'.txt'
;   skyapfile = 'skyaplist_'+dates+'.txt'
;   tracefile = 'tracelist_'+dates+'.txt'
;   crsplitfile = ['','crsplits_00mar30.txt','','','crsplits_00apr02.txt','']
;   crfile = 'crlist_'+dates+'.txt'
;   stdfile = 'stdlist_4.5_'+dates+'.txt'
;   stdfile2_5 = 'stdlist_2.5_'+dates+'.txt'
;   stdallfile = 'stdlist_4.5_00apr.txt'
;   stdallfile2_5 = 'stdlist_2.5_00apr.txt'
;   calibfile = 'caliblist_'+dates+'.txt'
;   tellfile = 'stdlist_'+dates+'.txt'
;   tellfits = 'telluric_'+dates+'.fits'
    
; ---------------------------------------------------------------------------
; pre-processing
; ---------------------------------------------------------------------------
    
;   iheader_check, root='a.'
;   ispec_makelists, ['a.28','a.29','a.30','a.31','a.32','a.33'], $
;     ['00mar29','00mar30','00mar31','00apr01','00apr02','00apr03'], /all, /overwrite, /gzip
    
;   skyflats = file_search('a.'+['280[4-5]','290[4-5]','300[4-5]','310[4-5]','320[5-6]','330[4-5]']+'.fits')
;   skyflatfile = 'skyflat_00apr.fits'
;   arm_flatcombine, skyflats, skyflatfile, [1208,5,1217,130], /display, $
;     psfile='qaplot_skyflat_00apr.ps'

;   imake_keywordfile, ['SCANLEN','POSANGLE','APERTURE'], comments=['scan length [arcsec]',$
;     'slit position angle [degrees]','slit aperture [arcsec]'], values=[0.0,90.0,2.5], $
;     root='a.01', headfile=headfile

;   headfile = 'header_keywords_00apr.dat'
;   iheader_keywords, headfile, /update, silent=silent

; ---------------------------------------------------------------------------
; carry out all the initial reductions
; ---------------------------------------------------------------------------
    
    iall, paramfile, procfile=procfile, doplot=doplot, tracefile=tracefile

; ---------------------------------------------------------------------------    
; combine cr-splits and cosmic-ray reject
; ---------------------------------------------------------------------------    

;   ibatch, paramfile[1], crsplitfile=crsplitfile[1], /find_pixshift, /debug, /crsplits

    iall, paramfile, crsplitfile=crsplitfile, crfile=crfile, /find_pixshift

; clean up residual cosmic rays

    ibatch, paramfile[2], crlist='ra.3011.fits', /crclean, objlim=0.5, iaxis=0
    ibatch, paramfile[2], crlist='ra.3021.fits', /crclean, objlim=0.5, iaxis=0
    ibatch, paramfile[3], crlist='ra.3127.fits', /crclean, objlim=0.5, iaxis=0
    ibatch, paramfile[4], crlist='ra.3221.fits', /crclean, objlim=0.4, iaxis=0
    ibatch, paramfile[5], crlist='ra.3334.fits', /crclean, objlim=0.5, iaxis=0

; repair LA_COSMIC pixels

    repair_la_cosmic_crpix, 'cra.3327.fits', 'repair.3327.dat', /repair

; ---------------------------------------------------------------------------    
; select sky apertures and sky subtract
; ---------------------------------------------------------------------------    

;   iskyselect, skyfile=skyfile, skyapfile=skyapfile, norder_sky=1, /skyremember

    iall, paramfile, skyapfile=skyapfile

; ---------------------------------------------------------------------------    
; make the sensitivity curves
; ---------------------------------------------------------------------------    

; generate sensitivity curves for each night (2.5")

    senslist = ['sens_2.5_00mar29.fits','sens_2.5_00mar30.fits','sens_2.5_00mar31.fits',$
      'sens_2.5_00apr01.fits','sens_2.5_00apr02.fits','sens_2.5_00apr03.fits']
    senstitle = '2000 '+['March 29','March 30','March 31','April 01','April 02','April 03']+' (2.5" Slit)'
    iall, paramfile, stdfile=stdfile2_5, senstitle=senstitle, sensinfo=sensinfo, $
      grey=2, slit_width=2.5, doplot=doplot, sensname=senslist
    
;   readcol, stdfile2_5[i], stdlist, format='A', /silent, comment='#'
;   info = isensfunc(i1dnames(stdlist,aperture=21),grey=2,sensname=senslist[i],$
;     senstitle=senstitle[i],doplot=doplot,/write,slit_width=2.5,zptshift=0.0)

; generate sensitivity curves for each night (4.5")

    senstitle = '2000 '+['March 29','March 30','March 31','April 01','April 02','April 03']
    iall, paramfile, stdfile=stdfile4_5, senstitle=senstitle, sensinfo=sensinfo, $
      grey=2, slit_width=4.5, doplot=doplot
    
;   senslist = ['sens_4.5_00mar29.fits','sens_4.5_00mar31.fits','sens_4.5_00apr01.fits',$
;     'sens_4.5_00apr02.fits','sens_4.5_00apr03.fits']
;   readcol, stdfile[i], stdlist, format='A', /silent, comment='#'
;   info = isensfunc(i1dnames(stdlist,aperture=21),grey=2,sensname=senslist[i],$
;     senstitle=senstitle[i],doplot=doplot,/write,slit_width=4.5,zptshift=0.0)

; make the combined sensitivity functions
    
    sensname = 'sens_2.5_00apr.fits' & senstitle='2000 April (2.5" Slit)'
    readcol, stdallfile2_5, stdlist, format='A', /silent, comment='#'
    ibatch, paramfile[0], stdlist=stdlist, sensname=sensname, /makesens, zptshift=0.0, $
      grey=2, sensinfo=sensinfo, senstitle=senstitle, slit_width=2.5, doplot=doplot, $
      wmapname='wmaplist_00apr.txt'

;   info = isensfunc(i1dnames(stdlist,aperture=21),grey=2,sensname=sensname,$
;     senstitle=senstitle,doplot=doplot,/write,slit_width=2.5,zptshift=0.0)

    sensname = 'sens_4.5_00apr.fits' & senstitle='2000 April (4.5" Slit)'
    readcol, stdallfile4_5, stdlist, format='A', /silent, comment='#'
    ibatch, paramfile[0], stdlist=stdlist, sensname=sensname, /makesens, zptshift=zptshift, $
      grey=2, sensinfo=sensinfo, senstitle=senstitle, slit_width=4.5, doplot=doplot, $
      wmapname='wmaplist_00apr.txt'

;   info = isensfunc(i1dnames(stdlist,aperture=21),grey=2,sensname=sensname,$
;     senstitle=senstitle,/doplot,/write,slit_width=4.5,zptshift=zptshift)

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

    sensname = 'sens_4.5_00apr.fits'
    iall, paramfile, fluxcalibfile=calibfile, calibfile=calibfile, $
      sensname=sensname;, tellfits=tellfits

; ---------------------------------------------------------------------------    
; miscellaneous tasks
; ---------------------------------------------------------------------------    

; compare the sensitivity curves

    label = [ ['29 March','Late clouds'],['30 March','Cloudy'],['31 March','Variable clouds'],$
      ['01 April','Clear'],['02 April','Clear'],['03 April','Clear'],['2.5" Slit',''] ]
    senslist = ['sens_4.5_00mar29.fits','sens_4.5_00mar30.fits','sens_4.5_00mar31.fits',$
      'sens_4.5_00apr01.fits','sens_4.5_00apr02.fits','sens_4.5_00apr03.fits','sens_2.5_00apr.fits']
    meansens = 'sens_4.5_00apr.fits'
    isenscompare, senslist, meansens=meansens, label=label, title='2000 April', $
      psname='qaplot_sens_compare_00apr.ps', /postscript

    label = [ ['29 March','Late clouds'],['30 March','Cloudy'],['31 March','Variable clouds'],$
      ['01 April','Clear'],['02 April','Clear'],['03 April','Clear'],['2.5" Slit',''],['4.5" Slit',''] ]
    loglist = ['qalog_sens_4.5_00mar29.log','qalog_sens_4.5_00mar30.log','qalog_sens_4.5_00mar31.log',$
      'qalog_sens_4.5_00apr01.log','qalog_sens_4.5_00apr02.log','qalog_sens_4.5_00apr03.log',$
      'qalog_sens_2.5_00apr.log','qalog_sens_4.5_00apr.log']
    iqaphoto, loglist, label=label, title='2000 April', psname='qaplot_photo_00apr.ps', /postscript
    
; ---------------------------------------------------------------------------    
; generate the web page
; ---------------------------------------------------------------------------    

    if keyword_set(web) then begin
       readcol, 'weblist_00apr.txt', weblist, format='A', comment='#', /silent
       ispec_webpage, '00apr', weblist=weblist, html_path=atlas_path(/dataweb), $
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

;   ibatch, paramfile[1], crlist='ra.3010.fits', /crclean, sigfrac=1.6, objlim=1.5 ; N3729
;   ibatch, paramfile[1], crlist='ra.3011.fits', /crclean, objlim=0.6 ; N4100
;   ibatch, paramfile[1], crlist='ra.3022.fits', /crclean, objlim=1.5, iaxis=1 ; N4010
;   ibatch, paramfile[2], crlist='ra.3123.fits', /crclean, objlim=1.2 ; U06917
;   ibatch, paramfile[2], crlist='ra.3125.fits', /crclean, objlim=1.2 ; U06816 
;   ibatch, paramfile[3], crlist='ra.3221.fits', /crclean, objlim=1.2, iaxis=1 ; U06969
;   ibatch, paramfile[3], crlist='ra.3227.fits', /crclean, sigfrac=1 ; N3893
;   ibatch, paramfile[4], crlist='ra.3314.fits', /crclean, objlim=1.5 ; IC749
;   ibatch, paramfile[4], crlist='ra.3334.fits', /crclean, objlim=1.5 ; N4111
    
; repair sky lines/stars that have been damaged by a cosmic rays

;   ibatch, paramfile[3], crlist='ra.3216_2.fits', /crclean, objlim=1.2 ; N3860 pos 6
;   rcube = rd2dspec('ra.3216_2.fits') ; i do not believe that it's a cosmic ray jm02nov21uofa
;   ccube = rd2dspec('cra.3216_2.fits')
;   ccube.image[697:724,40:45] = rcube.image[697:724,40:45]
;   ccube.mask[697:724,40:45] = fix(0)
;   wrt2dspec, ccube.fname, ccube.image, ccube.sigmamap, ccube.mask, *ccube.header
;   icleanup, rcube
;   icleanup, ccube

return
end
