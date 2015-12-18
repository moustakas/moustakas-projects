pro nov94_script, doplot=doplot, zptshift=zptshift, web=web

    paramfile = 'ibatch_94nov28.txt'
    procfile = 'objlist_94nov28.txt'
    skyfile = 'skylist_94nov28.txt'
    skyapfile = 'skyaplist_94nov28.txt'
    tracefile = 'tracelist_94nov28.txt'
    crsplitfile = 'crsplits_94nov28.txt'
    crfile = 'crlist_94nov28.txt'
    stdfile = 'stdlist_94nov28.txt'
    calibfile = 'caliblist_94nov28.txt'
    tellfile = 'stdlist_94nov28.txt'
    tellfits = 'telluric_94nov28.fits'
    
; ---------------------------------------------------------------------------
; pre-processing
; ---------------------------------------------------------------------------
    
;   iheader_check, root='n2?.'

;   skyflats = file_search('n28.'+['010[8-9]','011[0-3]']+'.fits')
;   skyflatfile = 'skyflat_94nov.fits'
;   arm_flatcombine, skyflats, skyflatfile, [1205,30,1217,285], /display, psfile='qaplot_skyflat_94nov.ps'
    
;   ispec_makelists, ['n27','n28'], ['94nov27','94nov28'], /all, /overwrite, $
;     gzip=0, gain=1.5, rdnoise=5.8, lampname='HeAr', extfile='kpnoextinct.dat', $
;     pscale=0.833333, minwave_guess=3600.0, minwave_out=3610.0, dwave=2.75, $
;     trim='3 1197 15 290', overscan='1205 1217 15 290', searchrad=3000.0, $
;     tracename='', badpixfile='badpix.dat'

; convert the dates to the new FITS standard

;   turner_new_fitsdate, /update

;   imake_keywordfile, ['SCANLEN','POSANGLE','APERTURE'], comments=['scan length [arcsec]',$
;     'slit position angle [degrees]','slit aperture [arcsec]'], values=[0.0,90.0,2.5], $
;     root='n2?.', headfile=headfile

;   headfile = 'header_keywords_94nov.dat'
;   iheader_keywords, headfile, /update, silent=silent

; ---------------------------------------------------------------------------
; carry out all the initial reductions
; ---------------------------------------------------------------------------
    
    iall, paramfile, procfile=procfile, doplot=doplot, tracefile=tracefile

; ---------------------------------------------------------------------------    
; combine cr-splits and cosmic-ray reject
; ---------------------------------------------------------------------------    

; REDUCE 0132/0132 (NGC942/3) both separately and together
    
;   ibatch, paramfile[0], crsplitfile=crsplitfile[0], /find_pixshift, /debug, /crsplits

    iall, paramfile, crsplitfile=crsplitfile, crfile=crfile, /find_pixshift

; ARP 230 drift -  the second exposure in this sequence (n28.0119) is
;                 too deviant from the other three images due to
;                 clouds and registration, so do not include it

    icrcombine, 'rn28.'+['0117','0120','0123']+'.fits', /find_pixshift, $
      user_pixshift=[0.0,-4.0,-2.0], crimname='rn28.0117_3.fits', /wfits;, /debug
    ibatch, paramfile[0], crlist='rn28.0117_3.fits', objlim=0.5, /crclean

; repair LA_COSMIC pixels

    repair_la_cosmic_crpix, 'crn28.0117_3.fits', 'repair.0117_3.dat', /repair
    
; ---------------------------------------------------------------------------    
; select sky apertures and sky subtract
; ---------------------------------------------------------------------------    
   
;   iskyselect, skyfile=skyfile, skyapfile=skyapfile, norder_sky=1, /skyremember
    iall, paramfile, skyapfile=skyapfile

; ---------------------------------------------------------------------------    
; make the sensitivity curve
; ---------------------------------------------------------------------------    

    senstitle = '1994 Nov 28'
    iall, paramfile, stdfile=stdfile, senstitle=senstitle, sensinfo=sensinfo, $
      slit_width=4.5, grey=2, doplot=doplot, searchrad=3000.0, zptshift=zptshift

;   senslist = 'sens_94nov28.fits'
;   readcol, stdfile, stdlist, format='A', /silent, comment='#'
;   info = isensfunc(i1dnames(stdlist,aperture=21.33),grey=2,sensname=senslist,$
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
       readcol, 'weblist_94nov.txt', weblist, format='A', comment='#', /silent
       ispec_webpage, '94nov', weblist=weblist, html_path=atlas_path(/dataweb), $
         html_only=html_only, /cleanpng
    endif

return 
end   
