pro may02_script, doplot=doplot, zptshift=zptshift, web=web
    
; RENAME all the files:
;
;   a.42??.fits --> a.45??.fits
;   a.43??.fits --> a.46??.fits
;   a.44??.fits --> a.47??.fits

    dateroot = '02may'
    dates = dateroot+['11','12','13']

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
;   ispec_makelists, ['a.45','a.46','a.47'], ['02may11','02may12','02may13'], /all, /overwrite, /gzip
    
;   skyflats = file_search('a.'+['460[3-5]','470[3-5]']+'.fits')
;   skyflatfile = 'skyflat_02may.fits'
;   arm_flatcombine, skyflats, skyflatfile, [1204,7,1217,138], /display, psfile='qaplot_skyflat_02may.ps'
    
;   imake_keywordfile, ['SCANLEN','POSANGLE','APERTURE'], comments=['scan length [arcsec]',$
;     'slit position angle [degrees]','slit aperture [arcsec]'], values=[0.0,90.0,2.5], $
;     root='a.01', headfile=headfile

;   headfile = 'header_keywords_02may.dat'
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

; NGC4559 West drift55
    icrcombine, 'ra.'+['4728','4729']+'.fits', /find_pixshift, $
      user_pixshift=[0.0,0.0], /wfits;, /debug
    ibatch, paramfile[2], crlist='ra.4728_2.fits', /crclean

; repair LA_COSMIC pixels

    repair_la_cosmic_crpix, 'cra.4513_2.fits', 'repair.4513_2.dat', /repair

; NGC4536 (4713) may need to be repaired near H-alpha, but since this
; nuclear spectrum is not in the atlas, forget about it
    
; ---------------------------------------------------------------------------    
; SINGS: stitch; select sky apertures and sky-subtract
; ---------------------------------------------------------------------------    

    sings2d_stitch, 'cra.'+['4508_2','4510_2']+'.fits', outname='wcra.4508_2_stitch.fits', $
      object='NGC4579 drift55', user_pixshift=[103.0], wmapname='wmap_ra.4503.idlsave', $
      refindx=0L, minwave=3640.0, dwave=2.75, scalefactor=[1.0,1.06], debug=1, wfits=1

    sings2d_stitch, 'cra.'+['4520_2','4522_2']+'.fits', outname='wcra.4520_2_stitch.fits', $
      object='NGC4826 drift55', user_pixshift=[101.0], wmapname='wmap_ra.4503.idlsave', $
      refindx=0L, minwave=3640.0, dwave=2.75, scalefactor=[1.0,1.08], debug=1, wfits=1

    sings2d_stitch, 'cra.'+['4618_3','4620_2']+'.fits', outname='wcra.4618_3_stitch.fits', $
      object='NGC4594 drift55', user_pixshift=[101.0], wmapname='wmap_ra.4615.idlsave', $
      refindx=0L, minwave=3640.0, dwave=2.75, scalefactor=[1.0,1.04], debug=1, wfits=1

; not a great profile match    
    sings2d_stitch, 'cra.'+['4626_2','4628_2']+'.fits', outname='wcra.4626_2_stitch.fits', $
      object='NGC4725 drift55', user_pixshift=[102.0], wmapname='wmap_ra.4615.idlsave', $
      refindx=0L, minwave=3640.0, dwave=2.75, scalefactor=[1.05,1.0], debug=1, wfits=1

    sings2d_stitch, 'cra.'+['4715_2','4717_2']+'.fits', outname='wcra.4715_2_stitch.fits', $
      object='NGC4536 drift55', user_pixshift=[102.0], wmapname='wmap_ra.4742.idlsave', $
      refindx=0L, minwave=3640.0, dwave=2.75, scalefactor=[1.0,1.02], debug=1, wfits=1

    sings2d_stitch, 'cra.'+['4726_2','4728_2']+'.fits', outname='wcra.4726_2_stitch.fits', $
      object='NGC4559 drift55', user_pixshift=[101.0], wmapname='wmap_ra.4742.idlsave', $
      refindx=0L, minwave=3640.0, dwave=2.75, scalefactor=[1.0,1.0], debug=1, wfits=1
    
;   readcol, 'skylist_stitch.txt', skylist, format='A', comment='#', /silent
;   iskyselect, skylist=skylist, skyapfile='skyaplist_stitch.txt'
    ibatch, paramfile[0], skyapfile='skyaplist_stitch.txt', skymethod=3L, /skysub
    
; ---------------------------------------------------------------------------    
; select sky apertures and sky subtract
; ---------------------------------------------------------------------------    

;   iskyselect, skyfile=skyfile, skyapfile=skyapfile, norder_sky=1, /skyremember
    iall, paramfile, skyapfile=skyapfile

; ---------------------------------------------------------------------------    
; make the sensitivity curves
; ---------------------------------------------------------------------------    

; generate sensitivity curves for each night (2.5")

    senstitle = '2002 '+['May 11','May 12','May 13']+' (2.5" Slit)'
    senslist = ['sens_2.5_02may11.fits','sens_2.5_02may12.fits','sens_2.5_02may13.fits']
    iall, paramfile, stdfile=stdfile2_5, senstitle=senstitle, sensinfo=sensinfo, $
      slit_width=2.5, grey=2, doplot=doplot, sensname=senslist

;   readcol, stdfile2_5[i], stdlist, format='A', /silent, comment='#'
;   info = isensfunc(i1dnames(stdlist,aperture=22.83),grey=2,sensname=senslist[i],$
;     senstitle=senstitle[i],doplot=doplot,/write,slit_width=2.5,zptshift=0.0)

; generate sensitivity curves for each night (4.5")

    senstitle = '2002 '+['May 11','May 12','May 13']
    iall, paramfile, stdfile=stdfile4_5, senstitle=senstitle, sensinfo=sensinfo, $
      slit_width=4.5, grey=2, doplot=doplot

;   senslist = ['sens_4.5_02may11.fits','sens_4.5_02may12.fits','sens_4.5_02may13.fits']
;   readcol, stdfile[i], stdlist, format='A', /silent, comment='#'
;   info = isensfunc(i1dnames(stdlist,aperture=22.83),grey=2,sensname=senslist[i],$
;     senstitle=senstitle[i],doplot=doplot,/write,slit_width=4.5,zptshift=0.0)

; make the combined sensitivity function
    
    sensname = 'sens_2.5_02may.fits' & senstitle='2002 May (2.5" Slit)'
    readcol, stdallfile2_5, stdlist, format='A', /silent, comment='#'
    ibatch, paramfile[0], stdlist=stdlist, sensname=sensname, /makesens, zptshift=0.0, $
      grey=2, sensinfo=sensinfo, senstitle=senstitle, slit_width=2.5, doplot=doplot, $
      wmapname='wmaplist_02may.txt'
    
;   info = isensfunc(i1dnames(stdlist,aperture=22.83),grey=2,sensname=sensname,$
;     senstitle=senstitle,doplot=doplot,/write,slit_width=2.5,zptshift=0.0)

    sensname = 'sens_4.5_02may.fits' & senstitle='2002 May (4.5" Slit)'
    readcol, stdallfile4_5, stdlist, format='A', /silent, comment='#'
    ibatch, paramfile[0], stdlist=stdlist, sensname=sensname, /makesens, zptshift=zptshift, $
      grey=2, sensinfo=sensinfo, senstitle=senstitle, slit_width=4.5, doplot=doplot, $
      wmapname='wmaplist_02may.txt'

;   info = isensfunc(i1dnames(stdlist,aperture=22.83),grey=2,sensname=sensname,$
;     senstitle=senstitle,doplot=doplot,/write,slit_width=4.5,zptshift=zptshift)

; ---------------------------------------------------------------------------    
; generate the telluric spectra
; ---------------------------------------------------------------------------    

    for i = 0L, n_elements(paramfile)-1L do begin

       readcol, tellfile[i], stdlist, format='A', /silent, comment='#'
       hotlist = i1dnames(stdlist,aperture=22)

       iconstruct_telluric, hotlist, tellmethod=2L, contmethod=1L, method=3L, $
         npoly=2, psname='qaplot_'+repstr(tellfits[i],'.fits','.ps'), tellfits=tellfits[i], $
         doplot=doplot, debug=0, /write

    endfor
    
; ---------------------------------------------------------------------------    
; flux-calibrate
; ---------------------------------------------------------------------------    
    
    sensname = 'sens_4.5_02may.fits'
    iall, paramfile, fluxcalibfile=calibfile, calibfile=calibfile, $
      sensname=sensname;, tellfits=tellfits

; the stitched spectra have been wavelength calibrated already, so
; just flux calibrate    
    
    readcol, 'skylist_stitch.txt', caliblist, format='A', comment='#', /silent
    ibatch, paramfile[0], caliblist='s'+caliblist, $ ; tellfits=tellfits[0], $
      sensname=sensname, /fluxcalibrate

; ---------------------------------------------------------------------------    
; miscellaneous tasks
; ---------------------------------------------------------------------------    

; compare the sensitivity curves

    label = [['11 May','Early cirrus'],['12 May','Cirrus'],['13 May','Clear'],['2.5" Slit',''] ]
    senslist = ['sens_4.5_02may11.fits','sens_4.5_02may12.fits','sens_4.5_02may13.fits','sens_2.5_02may.fits']
    meansens = 'sens_4.5_02may.fits'
    isenscompare, senslist, meansens=meansens, label=label, title='2002 May', $
      psname='qaplot_sens_compare_02may.ps', /postscript

    loglist = ['qalog_sens_4.5_02may11.log','qalog_sens_4.5_02may12.log',$
      'qalog_sens_4.5_02may13.log','qalog_sens_2.5_02may.log','qalog_sens_4.5_02may.log']
    label = [['11 May','Early cirrus'],['12 May','Cirrus'],['13 May','Clear'],$
      ['2.5" Slit',''],['4.5" slit','']]
    iqaphoto, loglist, label=label, title='2002 May', psname='qaplot_photo_02may.ps', /postscript

; ---------------------------------------------------------------------------    
; generate the web page
; ---------------------------------------------------------------------------    

    if keyword_set(web) then begin
       readcol, 'weblist_02may.txt', weblist, format='A', comment='#', /silent
       ispec_webpage, '02may', weblist=weblist, html_path=atlas_path(/dataweb), $
         html_only=html_only, /cleanpng
    endif

; ---------------------------------------------------------------------------    
; old code
; ---------------------------------------------------------------------------    

; fix individual frames for cosmic rays

;   ibatch, paramfile[0], crlist='ra.4513_2.fits', /crclean, sigfrac=1.0 ; N4569
;   ibatch, paramfile[0], crlist='ra.4526.fits', /crclean, objlim=1.8 ; N5033
;   ibatch, paramfile[1], crlist='ra.4620_2.fits', /crclean, objlim=1.1 ; N4594
;   ibatch, paramfile[1], crlist='ra.4628_2.fits', /crclean, objlim=1.3 ; N4725
;   ibatch, paramfile[2], crlist='ra.4715_2.fits', /crclean, objlim=1.6 ; N4536
;   ibatch, paramfile[2], crlist='ra.4722.fits', /crclean, objlim=1.5 ; TOL89
;   ibatch, paramfile[2], crlist='ra.4723.fits', /crclean, objlim=1.2 ; N4321
;   ibatch, paramfile[2], crlist='ra.4725.fits', /crclean, objlim=1.8 ; N4559
;   ibatch, paramfile[2], crlist='ra.4737.fits', /crclean, objlim=1.5, iaxis=1 ; N6946
;
;   rcube = rd2dspec('ra.4713.fits') ; i do not believe that it's a cosmic ray jm02oct22uofa
;   ccube = rd2dspec('cra.4713.fits')
;   ccube.image[1075:1083,61:69] = rcube.image[1075:1083,61:69]
;   ccube.mask[1075:1083,61:69] = fix(0)
;   wrt2dspec, ccube.fname, ccube.image, ccube.sigmamap, ccube.mask, *ccube.header
;   icleanup, rcube & icleanup, ccube
;
;   rcube = rd2dspec('ra.4724.fits') ; not a cosmic ray jm02dec3uofa
;   ccube = rd2dspec('cra.4724.fits')
;   ccube.image[1075:1081,66:73] = rcube.image[1075:1081,66:73]
;   ccube.mask[1075:1081,66:73] = fix(0)
;   wrt2dspec, ccube.fname, ccube.image, ccube.sigmamap, ccube.mask, *ccube.header
;   icleanup, rcube & icleanup, ccube

return
end    
