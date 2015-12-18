pro apr97_script, doplot=doplot, zptshift=zptshift, web=web

    paramfile = 'ibatch_97apr10.txt'
    procfile = 'objlist_97apr10.txt'
    skyfile = 'skylist_97apr10.txt'
    skyapfile = 'skyaplist_97apr10.txt'
    tracefile = 'tracelist_97apr10.txt'
    crsplitfile = 'crsplits_97apr10.txt'
    crfile = 'crlist_97apr10.txt'
    stdfile = 'stdlist_97apr10.txt'
    stdallfile = 'stdlist_94nov.txt'
    calibfile = 'caliblist_97apr10.txt'
    tellfile = 'stdlist_97apr10.txt'
    tellfits = 'telluric_97apr10.fits'
    
; ---------------------------------------------------------------------------
; pre-processing
; ---------------------------------------------------------------------------
    
;   iheader_check, root='a10.'
;
;   biasfile = 'a10.0101.fits'
;   domefile = 'a10.0102.fits'
;   skyflatfile = 'skyflat_97apr10.fits'
;
;   skyflats = file_search(['a10.010[6-9].fits','a10.0110.fits'])
;   arm_flatcombine, skyflats, skyflatfile, [1205,0,1217,269], /display
;   
;   ispec_makelists, 'a10', '97apr10', /all, /overwrite, gzip=0, gain=2.2, $
;     rdnoise=7.7, lampname='HeAr', extfile='kpnoextinct.dat', pscale=0.833333, $
;     minwave_guess=3560.0, minwave_out=3565.0, dwave=2.75, trim='5 1197 0 269', $
;     overscan='1205 1217 0 269', searchrad=3000.0, tracename='', $
;     badpixfile='badpix.dat', biasfile=biasfile, domefile=domefile, $
;     skyflatfile=skyflatfile

; convert the dates to the new FITS standard

;   turner_new_fitsdate, /update

;   imake_keywordfile, ['SCANLEN','POSANGLE','APERTURE'], comments=['scan length [arcsec]',$
;     'slit position angle [degrees]','slit aperture [arcsec]'], values=[0.0,90.0,2.5], $
;     root='a10.', headfile=headfile
 
;   headfile = 'header_keywords_97apr.dat'
    iheader_keywords, headfile, /update, silent=silent

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
    
; generate the sensitivity function

    senstitle = '1997 April 10'
    iall, paramfile, stdfile=stdfile, senstitle=senstitle, sensinfo=sensinfo, $
      slit_width=4.5, grey=2, doplot=doplot, searchrad=3000.0, zptshift=zptshift

;   readcol, stdfile, stdlist, format='A', /silent, comment='#'
;   senslist = 'sens_97apr10.fits'
;   info = isensfunc(i1dnames(stdlist,aperture=22.5),grey=2,sensname=senslist,$
;     senstitle=senstitle,/doplot,/write,slit_width=4.5,zptshift=zptshift)
    
; ---------------------------------------------------------------------------    
; generate a telluric spectrum
; ---------------------------------------------------------------------------    

    readcol, tellfile, stdlist, format='A', /silent, comment='#'
    hotlist = i1dnames(stdlist,aperture=22.5)

    iconstruct_telluric, hotlist, tellmethod=2L, contmethod=1L, method=3L, $
      npoly=2, psname='qaplot_'+repstr(tellfits,'.fits','.ps'), tellfits=tellfits, $
      doplot=doplot, debug=0, /write

; ---------------------------------------------------------------------------    
; flux-calibrate
; ---------------------------------------------------------------------------    
    
    iall, paramfile, fluxcalibfile=calibfile, calibfile=calibfile;, tellfits=tellfits

; ---------------------------------------------------------------------------    
; miscellaneous tasks
; ---------------------------------------------------------------------------    
        
; ---------------------------------------------------------------------------    
; generate the web page
; ---------------------------------------------------------------------------    

    if keyword_set(web) then begin
       readcol, 'weblist_97apr.txt', weblist, format='A', comment='#', /silent
       ispec_webpage, '97apr', weblist=weblist, html_path=atlas_path(/dataweb), $
         html_only=html_only, /cleanpng
    endif

return 
end   
