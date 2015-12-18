pro oct95_script, doplot=doplot, zptshift=zptshift, web=web

    dateroot = '95oct'
    dates = dateroot+['17','18']
    
    paramfile = 'ibatch_'+dates+'.txt'
    procfile = 'objlist_'+dates+'.txt'
    skyfile = 'skylist_'+dates+'.txt'
    skyapfile = 'skyaplist_'+dates+'.txt'
    tracefile = 'tracelist_'+dates+'.txt'
    crsplitfile = 'crsplits_'+dates+'.txt'
    crfile = 'crlist_'+dates+'.txt'
    stdfile = 'stdlist_'+dates+'.txt'
    stdallfile = 'stdlist_'+dateroot+'.txt'
    calibfile = 'caliblist_'+dates+'.txt'
    tellfile = 'stdlist_'+dates+'.txt'
    tellfits = 'telluric_'+dates+'.fits'
    
; ---------------------------------------------------------------------------
; pre-processing
; ---------------------------------------------------------------------------
    
;   iheader_check, root=['o17','o18']
 
;   biasfile1 = 'bias_95oct17.fits'
;   domefile1 = 'domeflat_95oct17.fits'
;   skyflatfile1 = 'skyflat_95oct17.fits'
;
;   biasfile2 = 'bias_95oct18.fits'
;   domefile2 = 'domeflat_95oct18.fits'
;   skyflatfile2 = 'skyflat_95oct18.fits'
;
;   biasfile = [biasfile1,biasfile2]
;   domefile = [domefile1,domefile2]
;   skyflatfile = [skyflatfile1,skyflatfile2]
    
; Oct 17
 
;   biases = [file_search('o17.010[2-9].fits'),file_search('o17.011[0-6].fits')]
;   arm_zerocombine, biases, biasfile1, /display, psfile='qaplot_bias_95oct17.ps'
;   
;   flats = [file_search('o17.011[7-9].fits'),file_search('o17.012[0-6].fits')]
;   arm_flatcombine, flats, domefile1, [1205,29,1217,272], /display, $
;     psfile='qaplot_domeflat_95oct17.ps'
;   
;   skyflats = [file_search('o17.0129.fits'),file_search('o17.013[0-3].fits')]
;   arm_flatcombine, skyflats, skyflatfile1, [1205,29,1217,272], /display, $
;     psfile='qaplot_skyflat_95oct17.ps'
 
; Oct 18
 
;   biases = [file_search('o18.010?.fits'),file_search('o18.0110.fits')]
;   arm_zerocombine, biases, biasfile2, /display, psfile='qaplot_bias_95oct18.ps'
;   
;   flats = [file_search('o18.011[1-9].fits'),file_search('o18.0120.fits')]
;   arm_flatcombine, flats, domefile2, [1205,29,1217,272], /display, $
;     psfile='qaplot_domeflat_95oct18.ps'
;   
;   skyflats = [file_search('o18.013[7-9].fits'),file_search('o18.014[0-2].fits')]
;   arm_flatcombine, skyflats, skyflatfile2, [1205,29,1217,272], /display, $
;     psfile='qaplot_skyflat_95oct18.ps'

;   ispec_makelists, ['o17','o18'], ['95oct17','95oct18'], /all, /overwrite, gzip=0, $
;     gain=2.2, rdnoise=7.7, lampname='HeAr', extfile='kpnoextinct.dat', pscale=0.833333, $
;     minwave_guess=3650.0, minwave_out=3650.0, dwave=2.75, trim='5 1195 29 272', $
;     overscan='1205 1217 29 272', searchrad=3000.0, tracename='', badpixfile='badpix.dat', $
;     biasfile=biasfile, domefile=domefile, skyflatfile=skyflatfile

; convert the dates to the new FITS standard

;   turner_new_fitsdate, /update
 
;   imake_keywordfile, ['SCANLEN','POSANGLE','APERTURE'], comments=['scan length [arcsec]',$
;     'slit position angle [degrees]','slit aperture [arcsec]'], values=[0.0,90.0,2.5], $
;     root=['o17','o18'], headfile=headfile
 
;   headfile = 'header_keywords_95oct.dat'
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

    ibatch, paramfile[0], crlist='ro17.0158_3.fits', /crclean, objlim=1.0, iaxis=0
    ibatch, paramfile[1], crlist='ro18.0146_3.fits', /crclean, objlim=1.0, iaxis=0

; ---------------------------------------------------------------------------    
; select sky apertures and sky subtract
; ---------------------------------------------------------------------------    
   
;   iskyselect, skyfile=skyfile, skyapfile=skyapfile, norder_sky=1, /skyremember
    iall, paramfile, skyapfile=skyapfile

; ---------------------------------------------------------------------------    
; make the sensitivity curves
; ---------------------------------------------------------------------------    
    
; generate the sensitivity function for each night (4.5")

    senstitle = '1995 '+['Oct 17','Oct 18']
    iall, paramfile, stdfile=stdfile, senstitle=senstitle, sensinfo=sensinfo, $
      slit_width=4.5, grey=2, doplot=doplot, searchrad=3000.0

;   readcol, stdfile[i], stdlist, format='A', /silent, comment='#'
;   senslist = 'sens_95oct'+['17','18']+'.fits'
;   info = isensfunc(i1dnames(stdlist,aperture=22.5),grey=2,sensname=senslist[i],$
;     senstitle=senstitle[i],/doplot,/write,slit_width=4.5,zptshift=0.0)

    sensname = 'sens_95oct.fits' & senstitle='1995 October (4.5" Slit)'
    readcol, stdallfile, stdlist, format='A', /silent, comment='#'
    ibatch, paramfile[0], stdlist=stdlist, sensname=sensname, /makesens, zptshift=zptshift, $
      grey=2, sensinfo=sensinfo, senstitle=senstitle, slit_width=4.5, doplot=doplot, $
      wmapname='wmaplist_95oct.txt', searchrad=3000.0

;   info = isensfunc(i1dnames(stdlist,aperture=20.33),grey=2,sensname=sensname,$
;     senstitle=senstitle,/doplot,/write,slit_width=4.5,zptshift=zptshift)
    
; ---------------------------------------------------------------------------    
; generate a telluric spectrum
; ---------------------------------------------------------------------------    

    for i = 0L, n_elements(paramfile)-1L do begin
       
       readcol, tellfile[i], stdlist, format='A', /silent, comment='#'
       hotlist = i1dnames(stdlist,aperture=20.33)

       iconstruct_telluric, hotlist, tellmethod=2L, contmethod=1L, method=3L, $
         npoly=2, psname='qaplot_'+repstr(tellfits[i],'.fits','.ps'), tellfits=tellfits[i], $
         doplot=doplot, debug=0, /write

    endfor
       
; ---------------------------------------------------------------------------    
; flux-calibrate
; ---------------------------------------------------------------------------    

    sensname = 'sens_95oct.fits'
    iall, paramfile, fluxcalibfile=calibfile, calibfile=calibfile, $
      sensname=sensname;, tellfits=tellfits

; ---------------------------------------------------------------------------    
; miscellaneous tasks
; ---------------------------------------------------------------------------    

; ---------------------------------------------------------------------------    
; generate the web page
; ---------------------------------------------------------------------------    

    if keyword_set(web) then begin
       readcol, 'weblist_95oct.txt', weblist, format='A', comment='#', /silent
       ispec_webpage, '95oct', weblist=weblist, html_path=atlas_path(/dataweb), $
         html_only=html_only, /cleanpng
    endif

return 
end   
