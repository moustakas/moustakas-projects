pro mar98_script, doplot=doplot, zptshift=zptshift, web=web

    dateroot = '98mar'
    dates = dateroot+['21']

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
    
;   iheader_check, root='a.01'

;   skyflats = file_search('a.010[5-7].fits')
;   arm_flatcombine, skyflats, 'skyflat_98mar21.fits', [1202,0,1217,125], /display, $
;     psfile='qaplot_skyflat_98mar21.ps'
 
;   ispec_makelists, 'a.01', '98mar21', /all, /overwrite, /gzip
    
;   imake_keywordfile, ['SCANLEN','POSANGLE','APERTURE'], comments=['scan length [arcsec]',$
;     'slit position angle [degrees]','slit aperture [arcsec]'], values=[0.0,90.0,2.5], $
;     root='a.01', headfile=headfile

;   headfile = 'header_keywords_98mar21.dat'
;   iheader_keywords, headfile, silent=silent, /update

; ---------------------------------------------------------------------------
; carry out all the initial reductions
; ---------------------------------------------------------------------------
    
    iall, paramfile, procfile=procfile, doplot=doplot, tracefile=tracefile

; ---------------------------------------------------------------------------    
; combine cr-splits and cosmic-ray reject
; ---------------------------------------------------------------------------    

;   ibatch, paramfile[0], crsplitfile=crsplitfile[0], /find_pixshift, /debug, /crsplits

    iall, paramfile, crsplitfile=crsplitfile, crfile=crfile, /find_pixshift
    
; ---------------------------------------------------------------------------    
; select sky apertures and sky subtract
; ---------------------------------------------------------------------------    
   
;   iskyselect, skyfile=skyfile, skyapfile=skyapfile, norder_sky=1, /skyremember

    iall, paramfile, skyapfile=skyapfile

; ---------------------------------------------------------------------------    
; make the sensitivity curves
; ---------------------------------------------------------------------------    

    sensname = 'sens_2.5_98mar21.fits' & senstitle = '1998 March (2.5" Slit)'
    readcol, stdfile2_5, stdlist, format='A', /silent, comment='#'
    ibatch, paramfile[0], stdlist=stdlist, sensname=sensname, /makesens, zptshift=0.0, $
      grey=2, sensinfo=sensinfo, senstitle=senstitle, slit_width=2.5, doplot=doplot

;   info = isensfunc(i1dnames(stdlist,aperture=21),grey=2,sensname=sensname,$
;     senstitle=senstitle,/doplot,/write,slit_width=2.5,zptshift=0.0)

    sensname = 'sens_4.5_98mar21.fits' & senstitle = '1998 March (4.5" Slit)'
    readcol, stdfile4_5, stdlist, format='A', /silent, comment='#'
    ibatch, paramfile[0], stdlist=stdlist, sensname=sensname, /makesens, zptshift=zptshift, $
      grey=2, sensinfo=sensinfo, senstitle=senstitle, slit_width=4.5, doplot=doplot

;   info = isensfunc(i1dnames(stdlist,aperture=21),grey=2,sensname=sensname,$
;     senstitle=senstitle,/doplot,/write,slit_width=4.5,zptshift=zptshift)

; ---------------------------------------------------------------------------    
; generate the telluric spectrum
; ---------------------------------------------------------------------------    

    readcol, tellfile, stdlist, format='A', /silent, comment='#'
    hotlist = i1dnames(stdlist,aperture=21)

    iconstruct_telluric, hotlist, tellmethod=2L, contmethod=1L, method=3L, $
      npoly=2, psname='qaplot_telluric_98mar21.ps', tellfits=tellfits, $
      doplot=doplot, debug=doplot, /write
    
; ---------------------------------------------------------------------------    
; flux-calibrate
; ---------------------------------------------------------------------------    
    
    sensname = 'sens_4.5_98mar21.fits'
    iall, paramfile, fluxcalibfile=calibfile, calibfile=calibfile, $
      sensname=sensname;, tellfits=tellfits

; ---------------------------------------------------------------------------    
; miscellaneous tasks
; ---------------------------------------------------------------------------    

; compare the sensitivity curves and compute the quality assurance on
; the photometry

    label = reform(['21 March','2.5" Slit',''],3,1)
    senslist = 'sens_2.5_98mar21.fits'
    meansens = 'sens_4.5_98mar21.fits'
    isenscompare, senslist, meansens=meansens, title='1998 March 21', $
      label=label, psname='qaplot_sens_compare_98mar.ps', /postscript
    
    loglist = ['qalog_sens_2.5_98mar21.log','qalog_sens_4.5_98mar21.log']
    label = [ ['21 March','2.5" Slit','Clear'],['21 March','4.5" Slit','Clear'] ]
    iqaphoto, loglist, label=label, title='1998 March 21', psname='qaplot_photo_98mar21.ps', /postscript

; ---------------------------------------------------------------------------    
; generate the web page
; ---------------------------------------------------------------------------    

    if keyword_set(web) then begin
       readcol, 'weblist_98mar.txt', weblist, format='A', comment='#', /silent
       ispec_webpage, '98mar', weblist=weblist, html_path=atlas_path(/dataweb), $
         html_only=html_only, /cleanpng
    endif

; ---------------------------------------------------------------------------    
; call the initial reductions using IBATCH
; ---------------------------------------------------------------------------    
;
;  ibatch, paramfile, /makeflat, checkillum=checkillum, doplot=doplot
;
;  readcol, procfile, proclist, format='A', /silent, comment='#'
;  ibatch, paramfile, proclist=proclist, /ccdproc, checkoverscan=checkoverscan
;
;  readcol, tracefile, tracelist, format='A', /silent, comment='#'
;  ibatch, paramfile, tracelist=tracelist, /distortion, doplot=doplot
;
;  ibatch, paramfile, /arcfit, debug=debug, doplot=doplot
;
;  ibatch, paramfile, crsplitfile=crsplitfile, /crsplits
;
;  ibatch, paramfile, skyapfile=skyapfile, /skysub, doplot=doplot
;  
;  readcol, crfile, crlist, format='A', /silent, comment='#'
;  ibatch, paramfile, crlist=crlist, /crclean

return
end
