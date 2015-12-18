pro apr96_script, doplot=doplot, zptshift=zptshift, web=web
    
    paramfile = 'ibatch_96apr22.txt'
    procfile = 'objlist_96apr22.txt'
    skyfile = 'skylist_96apr22.txt'
    skyapfile = 'skyaplist_96apr22.txt'
    tracefile = 'tracelist_96apr22.txt'
    crsplitfile = 'crsplits_96apr22.txt'
    crfile = 'crlist_96apr22.txt'
    stdfile = 'stdlist_96apr22.txt'
    stdallfile = 'stdlist_94nov.txt'
    calibfile = 'caliblist_96apr22.txt'
    tellfile = 'stdlist_96apr22.txt'
    tellfits = 'telluric_96apr22.fits'
    
; ---------------------------------------------------------------------------
; pre-processing
; ---------------------------------------------------------------------------
    
;   iheader_check, root='a22.'

; generate the average bias, dome flat, and sky flat

;   biasfile = 'bias_96apr22.fits'
;   domefile = 'domeflat_96apr22.fits'
;   skyflatfile = 'skyflat_96apr22.fits'

;   biases = file_search('a22.00[0-1]?.fits')
;   arm_zerocombine, biases, biasfile, /display, psfile='qaplot_bias_96apr22.ps'
    
;   flats = file_search('a22.002?.fits')
;   arm_flatcombine, flats, domefile, [1205,0,1217,269], /display, psfile='qaplot_domeflat_96apr22.ps'
    
;   skyflats = file_search('a22.010[1-4].fits')
;   arm_flatcombine, skyflats, skyflatfile, [1205,0,1217,269], /display, psfile='qaplot_skyflat_96apr22.ps'
    
;   ispec_makelists, 'a22', '96apr22', /all, /overwrite, gzip=0, gain=2.2, $
;     rdnoise=7.7, lampname='HeAr', extfile='kpnoextinct.dat', pscale=0.833333, $
;     minwave_guess=3560.0, minwave_out=3565.0, dwave=2.75, trim='5 1195 0 269', $
;     overscan='1205 1217 0 269', searchrad=3000.0, tracename='', $
;     badpixfile='badpix.dat', biasfile=biasfile, domefile=domefile, $
;     skyflatfile=skyflatfile

; convert the dates to the new FITS standard

;   turner_new_fitsdate, /update

;   imake_keywordfile, ['SCANLEN','POSANGLE','APERTURE'], comments=['scan length [arcsec]',$
;     'slit position angle [degrees]','slit aperture [arcsec]'], values=[0.0,90.0,2.5], $
;     root='a22.', headfile=headfile

;   headfile = 'header_keywords_96apr.dat'
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

; ---------------------------------------------------------------------------    
; select sky apertures and sky subtract
; ---------------------------------------------------------------------------    
   
;   iskyselect, skyfile=skyfile, skyapfile=skyapfile, norder_sky=1, /skyremember

    iall, paramfile, skyapfile=skyapfile

; ---------------------------------------------------------------------------    
; make the sensitivity curves
; ---------------------------------------------------------------------------    
    
; generate the sensitivity function

    senstitle = '1996 April 22'
    iall, paramfile, stdfile=stdfile, senstitle=senstitle, sensinfo=sensinfo, $
      slit_width=4.5, grey=2, doplot=doplot, searchrad=3000.0, zptshift=zptshift

;   readcol, stdfile, stdlist, format='A', /silent, comment='#'
;   senslist = 'sens_96apr22.fits'
;   info = isensfunc(i1dnames(stdlist,aperture=22.5),grey=2,sensname=senslist,$
;     senstitle=senstitle,/doplot,/write,slit_width=4.5,zptshift=zptshift)
    
; no telluric features are in the wavelength range, so do not generate
; a telluric absorption spectrum    

; ---------------------------------------------------------------------------    
; flux-calibrate
; ---------------------------------------------------------------------------    
    
    iall, paramfile, fluxcalibfile=calibfile, calibfile=calibfile

; ---------------------------------------------------------------------------    
; miscellaneous tasks
; ---------------------------------------------------------------------------    
        
; ---------------------------------------------------------------------------    
; generate the web page
; ---------------------------------------------------------------------------    

    if keyword_set(web) then begin
       readcol, 'weblist_96apr.txt', weblist, format='A', comment='#', /silent
       ispec_webpage, '96apr', weblist=weblist, html_path=atlas_path(/dataweb), $
         html_only=html_only, /cleanpng
    endif

return 
end   
