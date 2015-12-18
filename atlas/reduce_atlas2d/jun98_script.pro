pro jun98_script, doplot=doplot, zptshift=zptshift, web=web
; jm02jan3uofa

    dateroot = '98jun'
    dates = dateroot+['24','25','26']

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
;   ispec_makelists, ['a.08','a.09','a.10'], ['98jun24','98jun25','98jun26'], /all, /overwrite, /gzip

;   skyflatfile = 'skyflat_98jun.fits'
;   skyflats = file_search(['a.080[4-5].fits','a.090[4-5].fits','a.100[4-5].fits'])
;   arm_flatcombine, skyflats, skyflatfile, [1202,0,1217,128], /display, $
;     psfile='qaplot_skyflat_98jun.ps'
    
;   imake_keywordfile, ['SCANLEN','PA','APERTURE'], comments=['scan length [arcsec]',$
;     'slit position angle [degrees]','slit aperture [arcsec]'], values=[0.0,90.0,2.5], $
;     root='a.01', headfile=headfile

; correct some headers

;   list = ['a.0818.fits','a.0819.fits']
;   nlist = n_elements(list)
;   for i = 0L, nlist-1L do begin
;      h = headfits(list[i])
;      sxaddpar, h, 'IMAGETYP', 'dark'
;      sxaddpar, h, 'OBJECT', 'dark'
;      modfits, list[i], 0, h
;   endfor

;   headfile = 'header_keywords_98jun.dat'
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

; repair LA_COSMIC pixels

    repair_la_cosmic_crpix, ['cra.0928_2.fits','cra.0934_2.fits'], $
      ['repair.0928_2.dat','repair.0934_2.dat'], /repair
    
; ---------------------------------------------------------------------------    
; select sky apertures and sky subtract
; ---------------------------------------------------------------------------    
   
;   iskyselect, skyfile=skyfile, skyapfile=skyapfile, norder_sky=1, /skyremember

    iall, paramfile, skyapfile=skyapfile

; ---------------------------------------------------------------------------    
; make the sensitivity curves
; ---------------------------------------------------------------------------    

; generate sensitivity curves for each night (2.5")

    senslist = ['sens_4.5_98jun24.fits','sens_4.5_98jun25.fits','sens_4.5_98jun26.fits']
    senstitle = '1998 '+['June 24','June 25','June 26']+' (2.5" Slit)'
    iall, paramfile, stdfile=stdfile2_5, senstitle=senstitle, sensinfo=sensinfo, grey=2, $
      slit_width=2.5, doplot=doplot, sensname=senslist

;   readcol, stdfile2_5[i], stdlist, format='A', /silent, comment='#'
;   info = isensfunc(i1dnames(stdlist,aperture=21.5),grey=2,sensname=senslist[i],$
;     senstitle=senstitle[i],doplot=doplot,/write,slit_width=2.5,zptshift=0.0)

; generate sensitivity curves for each night (4.5")

    senstitle = '1998 '+['June 24','June 25','June 26']
    iall, paramfile, stdfile=stdfile4_5, senstitle=senstitle, sensinfo=sensinfo, $
      grey=2, slit_width=4.5, doplot=doplot

;   senslist = ['sens_4.5_98jun24.fits','sens_4.5_98jun25.fits','sens_4.5_98jun26.fits']
;   readcol, stdfile4_5[i], stdlist, format='A', /silent, comment='#'
;   info = isensfunc(i1dnames(stdlist,aperture=21.5),grey=2,sensname=senslist[i],$
;     senstitle=senstitle[i],doplot=doplot,/write,slit_width=4.5,zptshift=0.0)

; make the combined sensitivity function
    
    sensname = 'sens_2.5_98jun.fits' & senstitle='1998 June (2.5" Slit)'
    readcol, stdallfile2_5, stdlist, format='A', /silent, comment='#'
    ibatch, paramfile[0], stdlist=stdlist, sensname=sensname, /makesens, zptshift=0.0, $
      grey=2, sensinfo=sensinfo, senstitle=senstitle, slit_width=2.5, doplot=doplot, $
      wmapname='wmaplist_98jun.txt'

;   info = isensfunc(i1dnames(stdlist,aperture=21.5),grey=2,sensname=sensname,$
;     senstitle=senstitle,doplot=doplot,/write,slit_width=2.5,zptshift=0.0)

    sensname = 'sens_4.5_98jun.fits' & senstitle='1998 June (4.5" Slit)'
    readcol, stdallfile4_5, stdlist, format='A', /silent, comment='#'
    ibatch, paramfile[0], stdlist=stdlist, sensname=sensname, /makesens, zptshift=zptshift, $
      grey=2, sensinfo=sensinfo, senstitle=senstitle, slit_width=4.5, doplot=doplot, $
      wmapname='wmaplist_98jun.txt'

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

    sensname = 'sens_4.5_98jun.fits'
    iall, paramfile, fluxcalibfile=calibfile, calibfile=calibfile, $
      sensname=sensname;, tellfits=tellfits
    
; ---------------------------------------------------------------------------    
; miscellaneous tasks
; ---------------------------------------------------------------------------    

; compare the sensitivity curves

    label = [ ['24 June','Clear'],['25 June','Clear'],['26 June','Clear'],['2.5" Slit',''] ]
    senslist = ['sens_4.5_98jun24.fits','sens_4.5_98jun25.fits','sens_4.5_98jun26.fits','sens_2.5_98jun.fits']
    meansens = 'sens_4.5_98jun.fits'
    isenscompare, senslist, meansens=meansens, title='1998 June', $
      label=label, psname='qaplot_sens_compare_98jun.ps', /postscript

    loglist = ['qalog_sens_4.5_98jun24.log','qalog_sens_4.5_98jun25.log',$
      'qalog_sens_4.5_98jun26.log','qalog_sens_2.5_98jun.log','qalog_sens_4.5_98jun.log']
    label = [ ['24 June','Clear'],['25 June','Clear'],['26 June','Clear'],['2.5" Slit',''],['4.5" Slit',''] ]
    iqaphoto, loglist, label=label, title='1998 June', psname='qaplot_photo_98jun.ps', /postscript
    
; ---------------------------------------------------------------------------    
; generate the web page
; ---------------------------------------------------------------------------    

    if keyword_set(web) then begin
       readcol, 'weblist_98jun.txt', weblist, format='A', comment='#', /silent
       ispec_webpage, '98jun', weblist=weblist, html_path=atlas_path(/dataweb), $
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

;   ibatch, paramfile[0], crlist='ra.0833_2.fits', /crclean, objlim=3.0, sigclip=5.0 ; MRK307
;   ibatch, paramfile[1], crlist='ra.0932_2.fits', /crclean, objlim=1.8 ; MRK321
;   ibatch, paramfile[2], crlist='ra.1022_2.fits', /crclean, objlim=1.3 ; MRK809

return
end
