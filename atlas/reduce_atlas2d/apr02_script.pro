pro apr02_script, doplot=doplot, zptshift=zptshift, web=web

    dateroot = '02apr'
    dates = dateroot+['13','14','15']

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
;   ispec_makelists, ['a.42','a.43','a.44'], ['02apr13','02apr14','02apr15'], /all, /overwrite, /gzip

;   skyflats = file_search('a.'+['420[3-5]','430[3-5]']+'.fits')
;   skyflatfile = 'skyflat_02apr.fits'
;   arm_flatcombine, skyflats, skyflatfile, [1207,16,1217,151], /display, psfile='qaplot_skyflat_02apr.ps'
    
;   imake_keywordfile, ['SCANLEN','POSANGLE','APERTURE'], comments=['scan length [arcsec]',$
;     'slit position angle [degrees]','slit aperture [arcsec]'], values=[0.0,90.0,2.5], $
;     root='a.01', headfile=headfile

; update the headers    
    
;   headfile = 'header_keywords_02apr.dat'
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

; clean up residual cosmic rays

    ibatch, paramfile[0], crlist='ra.4215.fits', /crclean, sigclip=8.0, iaxis=0 ; LA_COSMIC is eating H-alpha
    ibatch, paramfile[0], crlist='ra.4264_2.fits', /crclean, sigclip=10.0, iaxis=0 ; LA_COSMIC is eating H-alpha
    ibatch, paramfile[2], crlist='ra.4415.fits', /crclean, objlim=1.0, iaxis=1

; repair LA_COSMIC pixels

    repair_la_cosmic_crpix, 'cra.4259_2.fits', 'repair.4259_2.dat', /repair

; ---------------------------------------------------------------------------    
; select sky apertures and sky subtract
; ---------------------------------------------------------------------------    

;   iskyselect, skyfile=skyfile, skyapfile=skyapfile, norder_sky=1, /skyremember

    iall, paramfile, skyapfile=skyapfile

; ---------------------------------------------------------------------------    
; make the sensitivity curves
; ---------------------------------------------------------------------------    

; generate sensitivity curves for each night (2.5")
    
    senstitle = '2002 '+['April 13','April 14','April 15']+'( 2.5" Slit)'
    senslist = ['sens_2.5_02apr13.fits','sens_2.5_02apr14.fits','sens_2.5_02apr15.fits']
    iall, paramfile, stdfile=stdfile2_5, senstitle=senstitle, sensinfo=sensinfo, $
      slit_width=2.5, grey=2, doplot=doplot, sensname=senslist

;   readcol, stdfile2_5[i], stdlist, format='A', /silent, comment='#'
;   info = isensfunc(i1dnames(stdlist,aperture=),grey=2,sensname=senslist[i],$
;     senstitle=senstitle[i],doplot=doplot,/write,slit_width=2.5,zptshift=0.0)

; generate sensitivity curves for each night (4.5")
    
    senstitle = '2002 '+['April 13','April 14','April 15']
    iall, paramfile, stdfile=stdfile4_5, senstitle=senstitle, sensinfo=sensinfo, $
      slit_width=4.5, grey=2, doplot=doplot

;   senslist = ['sens_4.5_02apr13.fits','sens_4.5_02apr14.fits','sens_4.5_02apr15.fits']
;   readcol, stdfile[i], stdlist, format='A', /silent, comment='#'
;   info = isensfunc(i1dnames(stdlist,aperture=),grey=2,sensname=senslist[i],$
;     senstitle=senstitle[i],doplot=doplot,/write,slit_width=4.5,zptshift=0.0)

; make the combined sensitivity function
    
    sensname = 'sens_2.5_02apr.fits' & senstitle='2002 April (2.5" Slit)'
    readcol, stdallfile2_5, stdlist, format='A', /silent, comment='#'
    ibatch, paramfile[0], stdlist=stdlist, sensname=sensname, /makesens, zptshift=0.0, $
      grey=2, sensinfo=sensinfo, senstitle=senstitle, slit_width=2.5, doplot=doplot, $
      wmapname='wmaplist_02apr.txt'

;   info = isensfunc(i1dnames(stdlist,aperture=22.67),grey=2,sensname=sensname,$
;     senstitle=senstitle,doplot=doplot,/write,slit_width=2.5,zptshift=0.0)

    sensname = 'sens_4.5_02apr.fits' & senstitle='2002 April (4.5" Slit)'
    readcol, stdallfile4_5, stdlist, format='A', /silent, comment='#'
    ibatch, paramfile[0], stdlist=stdlist, sensname=sensname, /makesens, zptshift=zptshift, $
      grey=2, sensinfo=sensinfo, senstitle=senstitle, slit_width=4.5, doplot=doplot, $
      wmapname='wmaplist_02apr.txt'

;   info = isensfunc(i1dnames(stdlist,aperture=22.67),grey=2,sensname=sensname,$
;     senstitle=senstitle,doplot=doplot,/write,slit_width=4.5,zptshift=zptshift)

; ---------------------------------------------------------------------------    
; generate the telluric spectra
; ---------------------------------------------------------------------------    

    for i = 0L, n_elements(paramfile)-1L do begin

       readcol, tellfile[i], stdlist, format='A', /silent, comment='#'
       hotlist = i1dnames(stdlist,aperture=22.67)

       iconstruct_telluric, hotlist, tellmethod=2L, contmethod=1L, method=3L, $
         npoly=2, psname='qaplot_'+repstr(tellfits[i],'.fits','.ps'), tellfits=tellfits[i], $
         doplot=doplot, debug=0, /write

    endfor
    
; ---------------------------------------------------------------------------    
; flux-calibrate
; ---------------------------------------------------------------------------    
    
    sensname = 'sens_4.5_02apr.fits'
    iall, paramfile, fluxcalibfile=calibfile, calibfile=calibfile, $
      sensname=sensname;, tellfits=tellfits

; ---------------------------------------------------------------------------    
; miscellaneous tasks
; ---------------------------------------------------------------------------    

; compare the sensitivity curves

    label = [['13 April','Early clouds'],['14 April','Clear, windy'],['15 April','Cloudy, windy'],['2.5" Slit',''] ]
    senslist = ['sens_4.5_02apr13.fits','sens_4.5_02apr14.fits','sens_4.5_02apr15.fits','sens_2.5_02apr.fits']
    meansens = 'sens_4.5_02apr.fits'
    isenscompare, senslist, meansens=meansens, label=label, psname='qaplot_sens_compare_02apr.ps', $
      title='2002 April', /postscript

    loglist = ['qalog_sens_4.5_02apr13.log','qalog_sens_4.5_02apr14.log',$
      'qalog_sens_4.5_02apr15.log','qalog_sens_2.5_02apr.log','qalog_sens_4.5_02apr.log']
    label = [['13 April','Early clouds'],['14 April','Clear, windy'],$
      ['15 April','Cloudy, windy'],['2.5" Slit',''],['4.5" slit','']]
    iqaphoto, loglist, label=label, title='2002 April', psname='qaplot_photo_02apr.ps', /postscript
    
; ---------------------------------------------------------------------------    
; generate the web page
; ---------------------------------------------------------------------------    

    if keyword_set(web) then begin
       readcol, 'weblist_02apr.txt', weblist, format='A', comment='#', /silent
       ispec_webpage, '02apr', weblist=weblist, html_path=atlas_path(/dataweb), $
         html_only=html_only, /cleanpng
    endif

; ---------------------------------------------------------------------------    
; old code
; ---------------------------------------------------------------------------    
    
; fix specific images for cosmic rays

;    ibatch, paramfile[0], crlist='ra.4253.fits', /crclean, sigfrac=1.0, objlim=1.5 ; U09081
;    ibatch, paramfile[0], crlist='ra.4259_2.fits', /crclean, sigfrac=1.0 ; MRK848
;    
;    ibatch, paramfile[0], crlist='ra.4215.fits', /crclean, objlim=1.5 ; N3310
;    rcube = rd2dspec('ra.4215.fits') ; not cosmic rays jm02nov21uofa
;    ccube = rd2dspec('cra.4215.fits')
;    ccube.image[1078:1090,54:62] = rcube.image[1078:1090,54:62] & ccube.mask[1078:1090,54:62] = fix(0)
;    ccube.image[502:506,57:64] = rcube.image[502:506,57:64] & ccube.mask[502:506,57:64] = fix(0)
;    ccube.image[515:532,65:67] = rcube.image[515:532,65:67] & ccube.mask[515:532,65:67] = fix(0)
;    wrt2dspec, ccube.fname, ccube.image, ccube.sigmamap, ccube.mask, *ccube.header
;    icleanup, rcube & icleanup, ccube
;
;    ibatch, paramfile[0], crlist='ra.4242.fits', /crclean, objlim=1.5
;    rcube = rd2dspec('ra.4242.fits') ; not a cosmic ray jm02nov21uofa
;    ccube = rd2dspec('cra.4242.fits')
;    ccube.image[249:259,65:67] = rcube.image[249:259,65:67] & ccube.mask[249:259,65:67] = fix(0)
;    wrt2dspec, ccube.fname, ccube.image, ccube.sigmamap, ccube.mask, *ccube.header
;    icleanup, rcube & icleanup, ccube
;
;    ibatch, paramfile[0], crlist='ra.4264_2.fits', /crclean, objlim=2.5, sigfrac=1.2 ; N6052
;    rcube = rd2dspec('ra.4264_2.fits') ; not a cosmic ray jm02nov21uofa
;    ccube = rd2dspec('cra.4264_2.fits')
;    ccube.image[1114:1116,63:64] = rcube.image[1114:1116,63:64] & ccube.mask[1114:1116,63:64] = fix(0)
;    wrt2dspec, ccube.fname, ccube.image, ccube.sigmamap, ccube.mask, *ccube.header
;    icleanup, rcube & icleanup, ccube
;
;    ibatch, paramfile[1], crlist='ra.4327.fits', /crclean, objlim=1.5 ; N3953
;    ibatch, paramfile[1], crlist='ra.4356.fits', /crclean, objlim=2.2 ; HD167451
;    ibatch, paramfile[2], crlist='ra.4415.fits', /crclean, objlim=1.0, sigclip=4.0 ; N3628
    
return
end
