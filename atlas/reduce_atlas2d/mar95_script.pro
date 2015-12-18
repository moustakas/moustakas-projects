pro mar95_script, doplot=doplot, zptshift=zptshift, web=web

    dateroot = '95mar'
    dates = dateroot+['26','27']
    
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
    
;   iheader_check, root=['m26','m27']

;   biasfile1 = 'm26.Zero.fits'
;   domefile1 = 'm26.Flat.fits'
;   skyflatfile1 = 'skyflat_95mar26.fits'
;
;   biasfile2 = 'm27.Zero.fits'
;   domefile2 = 'm27.Dome.fits'
;   skyflatfile2 = 'skyflat_95mar27.fits'
;
;   biasfile = [biasfile1,biasfile2]
;   domefile = [domefile1,domefile2]
;   skyflatfile = [skyflatfile1,skyflatfile2]
    
; Mar 26    

;   skyflats = file_search('m26.010[1-7].fits')
;   arm_flatcombine, skyflats, skyflatfile1, [1205,0,1217,269], /display

; Mar 27
    
;   skyflats = file_search('m27.010[3-6].fits')
;   arm_flatcombine, skyflats, skyflatfile2, [1205,0,1217,269], /display

;   ispec_makelists, ['m26','m27'], ['95mar26','95mar27'], /all, /overwrite, gzip=0, $
;     gain=2.2, rdnoise=7.7, lampname='HeAr', extfile='kpnoextinct.dat', pscale=0.833333, $
;     minwave_guess=3560.0, minwave_out=3565.0, dwave=2.75, trim='5 1195 0 269', $
;     overscan='1205 1217 0 269', searchrad=3000.0, tracename='', badpixfile='badpix.dat', $
;     biasfile=biasfile, domefile=domefile, skyflatfile=skyflatfile

; convert the dates to the new FITS standard

;   turner_new_fitsdate, /update
 
;   imake_keywordfile, ['SCANLEN','POSANGLE','APERTURE'], comments=['scan length [arcsec]',$
;     'slit position angle [degrees]','slit aperture [arcsec]'], values=[0.0,90.0,2.5], $
;     root=['m26','m27'], headfile=headfile
 
;   headfile = 'header_keywords_95mar.dat'
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
    
; generate the sensitivity function for each night (4.5")

    senstitle = '1995 '+['Mar 26','Mar 27']
    iall, paramfile, stdfile=stdfile, senstitle=senstitle, sensinfo=sensinfo, $
      slit_width=4.5, grey=2, doplot=doplot, searchrad=3000.0

;   readcol, stdfile[i], stdlist, format='A', /silent, comment='#'
;   senslist = 'sens_95mar'+['26','27']+'.fits'
;   info = isensfunc(i1dnames(stdlist,aperture=22.5),grey=2,sensname=senslist[i],$
;     senstitle=senstitle[i],/doplot,/write,slit_width=4.5,zptshift=0.0)

    sensname = 'sens_95mar.fits' & senstitle='1995 March (4.5" Slit)'
    readcol, stdallfile, stdlist, format='A', /silent, comment='#'
    ibatch, paramfile[0], stdlist=stdlist, sensname=sensname, /makesens, zptshift=zptshift, $
      grey=2, sensinfo=sensinfo, senstitle=senstitle, slit_width=4.5, doplot=doplot, $
      wmapname='wmaplist_95mar.txt', searchrad=3000.0

;   info = isensfunc(i1dnames(stdlist,aperture=22.5),grey=2,sensname=sensname,$
;     senstitle=senstitle,/doplot,/write,slit_width=4.5,zptshift=zptshift)
    
; ---------------------------------------------------------------------------    
; generate a telluric spectrum
; ---------------------------------------------------------------------------    

    for i = 0L, n_elements(paramfile)-1L do begin
       
       readcol, tellfile[i], stdlist, format='A', /silent, comment='#'
       hotlist = i1dnames(stdlist,aperture=22.5)

       iconstruct_telluric, hotlist, tellmethod=2L, contmethod=1L, method=3L, $
         npoly=2, psname='qaplot_'+repstr(tellfits[i],'.fits','.ps'), tellfits=tellfits[i], $
         doplot=doplot, debug=0, /write

    endfor
       
; ---------------------------------------------------------------------------    
; flux-calibrate
; ---------------------------------------------------------------------------    
    
    sensname = 'sens_95mar.fits'
    iall, paramfile, fluxcalibfile=calibfile, calibfile=calibfile, $
      sensname=sensname;, tellfits=tellfits

; ---------------------------------------------------------------------------    
; miscellaneous tasks
; ---------------------------------------------------------------------------    

; ---------------------------------------------------------------------------    
; generate the web page
; ---------------------------------------------------------------------------    

    if keyword_set(web) then begin
       readcol, 'weblist_95mar.txt', weblist, format='A', comment='#', /silent
       ispec_webpage, '95mar', weblist=weblist, html_path=atlas_path(/dataweb), $
         html_only=html_only, /cleanpng
    endif

return 
end   
