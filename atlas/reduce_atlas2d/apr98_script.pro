pro apr98_script, doplot=doplot, zptshift=zptshift, web=web
; jm02jan3uofa

    dateroot_apr = '98apr'
    dateroot_may = '98may'
    dates = [dateroot_apr+['27','28','29','30'],dateroot_may+['01','02']]

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
    
; ---------------------------------------------------------------------------
; pre-processing
; ---------------------------------------------------------------------------
    
;   iheader_check, root='a.'
;   ispec_makelists, ['a.02','a.03','a.04','a.05','a.06','a.07'], $
;     ['98apr27','98apr28','98apr29','98apr30','98may01','98may02'], /all, /overwrite, /gzip
    
;   skyflatfile = 'skyflat_98apr.fits'
;   skyflats = file_search('a.'+['023[6-7]','030[5-7]','040[5-6]','050[7-8]','060[6-7]','070[6-7]']+'.fits')
;   arm_flatcombine, skyflats, skyflatfile, [1202,0,1217,128], /display, $
;     psfile='qaplot_skyflat_98apr.ps'
    
;   imake_keywordfile, ['SCANLEN','PA','APERTURE'], comments=['scan length [arcsec]',$
;     'slit position angle [degrees]','slit aperture [arcsec]'], values=[0.0,90.0,2.5], $
;     root='a.01', headfile=headfile

; correct some headers

;   list = ['a.0201.fits','a.0202.fits','a.0235.fits','a.0303.fits','a.0304.fits','a.0344.fits']
;   nlist = n_elements(list)
;   for i = 0L, nlist-1L do begin
;      h = headfits(list[i])
;      sxaddpar, h, 'IMAGETYP', 'comp'
;      modfits, list[i], 0, h
;   endfor
;
;   h = headfits('a.0224.fits')
;   sxaddpar, h, 'OBJECT', 'HZ 44 4.5'
;   modfits, 'a.0224.fits', 0, h
    
; fix one image

;   rowfix_0415, /wfits
    
; update the headers    
    
;   headfile = 'header_keywords_98apr.dat'
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

; MRK5 drift
    icrcombine, 'ra.'+['0512','0513']+'.fits', /find_pixshift, $
      user_pixshift=[0.0,1.0], /wfits;, /debug
    ibatch, paramfile[3], crlist='ra.0512_2.fits', /crclean

; NGC3198 drift
    icrcombine, 'ra.'+['0526','0527']+'.fits', /find_pixshift, $
      user_pixshift=[0.0,-5.0], /wfits;, /debug
    ibatch, paramfile[3], crlist='ra.0526_2.fits', /crclean

; NGC4713 drift
    icrcombine, 'ra.'+['0625','0626']+'.fits', /find_pixshift, $
      user_pixshift=[0.0,1.0], /wfits;, /debug
    ibatch, paramfile[4], crlist='ra.0625_2.fits', /crclean

; LA_COSMIC is over-zealous on this object 

    ibatch, paramfile[0], crlist='ra.0213_2.fits', /crclean, sigclip=10, iaxis=0

; ---------------------------------------------------------------------------    
; select sky apertures and sky subtract
; ---------------------------------------------------------------------------    

;   iskyselect, skyfile=skyfile, skyapfile=skyapfile, norder_sky=1, /skyremember

    iall, paramfile, skyapfile=skyapfile

; ---------------------------------------------------------------------------    
; make the sensitivity curves
; ---------------------------------------------------------------------------    

; generate sensitivity curves for each night (2.5")

    senstitle = '1998 '+['April 27','April 28','April 29','April 30','May 01','May 02']+'( 2.5" Slit)'
    senslist = ['sens_2.5_98apr27.fits','sens_2.5_98apr28.fits','sens_2.5_98apr29.fits',$
      'sens_2.5_98apr30.fits','sens_2.5_98may01.fits','sens_2.5_98may02.fits']
    iall, paramfile, stdfile=stdfile2_5, senstitle=senstitle, sensinfo=sensinfo, grey=2, $
      slit_width=2.5, doplot=doplot, sensname=senslist

;   readcol, stdfile2_5[i], stdlist, format='A', /silent, comment='#'
;   info = isensfunc(i1dnames(stdlist,aperture=21.5),grey=2,sensname=senslist[i],$
;     senstitle=senstitle[i],doplot=doplot,/write,slit_width=4.5,zptshift=0.0)

; generate sensitivity curves for each night (4.5")

    senstitle = '1998 '+['April 27','April 28','April 29','April 30','May 01','May 02']
    iall, paramfile, stdfile=stdfile4_5, senstitle=senstitle, sensinfo=sensinfo, grey=2, $
      slit_width=4.5, doplot=doplot

;   senslist = ['sens_4.5_98apr27.fits','sens_4.5_98apr28.fits','sens_4.5_98apr29.fits',$
;     'sens_4.5_98apr30.fits','sens_4.5_98may01.fits','sens_4.5_98may02.fits']
;   readcol, stdfile[i], stdlist, format='A', /silent, comment='#'
;   info = isensfunc(i1dnames(stdlist,aperture=21.5),grey=2,sensname=senslist[i],$
;     senstitle=senstitle[i],doplot=doplot,/write,slit_width=4.5,zptshift=0.0)

; drift-scanned standard star.  unknown scan length!
    
;;  ibatch, paramfile[3], stdlist='cra.0541.fits', sensname='sens_drift_98apr.fits', /makesens, $
;;    grey=2, sensinfo=sensinfo, senstitle=senstitle, slit_width=20.0, doplot=doplot
    
; make the combined sensitivity function
    
    sensname = 'sens_2.5_98apr.fits' & senstitle = '1998 April (2.5" Slit)'
    readcol, stdallfile2_5, stdlist, format='A', /silent, comment='#'
    ibatch, paramfile[0], stdlist=stdlist, sensname=sensname, /makesens, zptshift=0.0, $
      grey=2, sensinfo=sensinfo, senstitle=senstitle, slit_width=2.5, doplot=doplot, $
      wmapname='wmaplist_98apr.txt'

;   info = isensfunc(i1dnames(stdlist,aperture=21.5),grey=2,sensname=sensname,$
;     senstitle=senstitle,doplot=doplot,/write,slit_width=2.5,zptshift=0.0)

    sensname = 'sens_4.5_98apr.fits' & senstitle = '1998 April (4.5" Slit)'
    readcol, stdallfile4_5, stdlist, format='A', /silent, comment='#'
    ibatch, paramfile[0], stdlist=stdlist, sensname=sensname, /makesens, zptshift=zptshift, $
      grey=2, sensinfo=sensinfo, senstitle=senstitle, slit_width=4.5, doplot=doplot, $
      wmapname='wmaplist_98apr.txt'

;   info = isensfunc(i1dnames(stdlist,aperture=21.5),grey=2,sensname=sensname,$
;     senstitle=senstitle,doplot=doplot,/write,slit_width=4.5,zptshift=zptshift)

; ---------------------------------------------------------------------------    
; generate the telluric spectra
; ---------------------------------------------------------------------------    

    for i = 0L, n_elements(paramfile)-1L do begin

       readcol, tellfile[i], stdlist, format='A', /silent, comment='#'
       hotlist = i1dnames(stdlist,aperture=21.5)

       iconstruct_telluric, hotlist, tellmethod=2L, contmethod=1L, method=3L, $
         npoly=2, psname='qaplot_'+repstr(tellfits[i],'.fits','.ps'), tellfits=tellfits[i], $
         doplot=doplot, debug=0, /write

    endfor

; ---------------------------------------------------------------------------    
; flux-calibrate
; ---------------------------------------------------------------------------    

    sensname = 'sens_4.5_98apr.fits'
    iall, paramfile, fluxcalibfile=calibfile, calibfile=calibfile, sensname=sensname;, tellfits=tellfits

; ---------------------------------------------------------------------------    
; miscellaneous tasks
; ---------------------------------------------------------------------------    

; compare the sensitivity curves

    label = [ ['27 April','Clear'],['28 April','Clear'],['29 April','Early clouds'],$
      ['30 April','Clear'],['01 May','Late thin clouds'],['02 May','Clear'],['2.5" Slit',''] ]
    senslist = ['sens_4.5_98apr27.fits','sens_4.5_98apr28.fits','sens_4.5_98apr29.fits',$
      'sens_4.5_98apr30.fits','sens_4.5_98may01.fits','sens_4.5_98may02.fits','sens_2.5_98apr.fits']
    meansens = 'sens_4.5_98apr.fits'
    isenscompare, senslist, meansens=meansens, title='1998 April', $
      label=label, psname='qaplot_sens_compare_98apr.ps', /postscript

    loglist = ['qalog_sens_4.5_98apr27.log','qalog_sens_4.5_98apr28.log',$
      'qalog_sens_4.5_98apr29.log','qalog_sens_4.5_98apr30.log','qalog_sens_4.5_98may01.log',$
      'qalog_sens_4.5_98may02.log','qalog_sens_2.5_98apr.log','qalog_sens_4.5_98apr.log']
    label = [ ['27 April','Clear'],['28 April','Clear?'],['29 April','Cloudy?'],$
      ['30 April','Clear'],['01 May','Late thin clouds'],['02 May','Clear'],['2.5" Slit',''],['4.5" Slit',''] ]
    iqaphoto, loglist, label=label, title='1998 April', psname='qaplot_photo_98apr.ps', /postscript
    
; ---------------------------------------------------------------------------    
; generate the web page
; ---------------------------------------------------------------------------    

    if keyword_set(web) then begin
       readcol, 'weblist_98apr.txt', weblist, format='A', comment='#', /silent
       ispec_webpage, '98apr', weblist=weblist, html_path=atlas_path(/dataweb), $
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
    
; repair specific images for cosmic rays

;   ibatch, paramfile[0], crlist='ra.0213_2.fits', /crclean, objlim=1.5, sigclip=6.0, sigfrac=1.2 ; N5236 center
;   ibatch, paramfile[1], crlist='ra.0319_2.fits', /crclean, sigfrac=1.2, objlim=1.8 ; N3351 center
;   ibatch, paramfile[3], crlist='ra.0514_2.fits', /crclean, objlim=1.0, sigclip=4.0 ; U04483
;   ibatch, paramfile[4], crlist='ra.0625_2.fits', /crclean, objlim=1.5 ; N4713
;   ibatch, paramfile[4], crlist='ra.0639.fits', /crclean, sigfrac=1.4, objlim=1.8 ; MRK496AB
    
return
end
