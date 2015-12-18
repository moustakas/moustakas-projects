pro nov01_script, doplot=doplot, zptshift=zptshift, web=web

    dateroot = '01nov'
    dates = dateroot+['10','11','12','13']

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
;   ispec_makelists, ['a.34','a.35','a.36','a.37'], ['01nov10','01nov11','01nov12','01nov13'], $
;      /all, /overwrite, /gzip
    
;   skyflats = file_search('a.'+['340[4-6]','350[1-3]']+'.fits')
;   skyflatfile = 'skyflat_01nov.fits'
;   arm_flatcombine, skyflats, skyflatfile, [1206,3,1218,132], /display, $
;     psfile='qaplot_skyflat_01nov.ps'

;   imake_keywordfile, ['SCANLEN','POSANGLE','APERTURE'], comments=['scan length [arcsec]',$
;     'slit position angle [degrees]','slit aperture [arcsec]'], values=[0.0,90.0,2.5], $
;     root='a.01', headfile=headfile

; correct some headers

;   h = headfits('a.3516.fits')
;   sxaddpar, h, 'OBJECT', 'N7678 nuc'
;   modfits, 'a.3516.fits', 0, h
;   
;   h = headfits('a.3533.fits')
;   sxaddpar, h, 'OBJECT', 'N1482 east knot'
;   modfits, 'a.3533.fits', 0, h

; update the headers    
    
;   headfile = 'header_keywords_01nov.dat'
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

    repair_la_cosmic_crpix, 'cra.3516.fits', 'repair.3516.dat', /repair

; ---------------------------------------------------------------------------    
; SINGS: NGC0628 (3418-3421) and NGC0925 (3518-3521) must be stitched;
; select sky apertures and sky-subtract these spectra separately
; ---------------------------------------------------------------------------    

    sings2d_stitch, 'cra.'+['3418_2','3420_2']+'.fits', outname='wcra.3418_2_stitch.fits', $
      object='NGC0628 drift55', user_pixshift=[109.0], wmapname='wmap_ra.3447.idlsave', $
      refindx=0L, minwave=3630.0, dwave=2.75, scalefactor=[1.0,1.02], debug=1, wfits=1

; totally wrong stitch! remove from further analysis
;   sings2d_stitch, 'cra.'+['3438_2','3440_2']+'.fits', outname='wcra.3438_2_stitch.fits', $
;     object='HoII drift55', user_pixshift=[108.0], wmapname='wmap_ra.3447.idlsave', $
;     refindx=0L, minwave=3630.0, dwave=2.75, scalefactor=[1.0,1.0], stitchrange=[1.5,3]*1E5, debug=1, wfits=0

    sings2d_stitch, 'cra.'+['3518_2','3520_2']+'.fits', outname='wcra.3518_2_stitch.fits', $
      object='NGC0925 drift55', user_pixshift=[102.0], wmapname='wmap_ra.3508.idlsave', $
      refindx=0L, minwave=3630.0, dwave=2.75, scalefactor=[1.01,1.0], debug=1, wfits=1

;   readcol, 'skylist_stitch.txt', skylist, format='A', comment='#', /silent
;   iskyselect, skylist=skylist, skyapfile='skyaplist_stitch.txt'
    ibatch, paramfile[0], skyapfile='skyaplist_stitch.txt', skymethod=3L, /skysub
    
; ---------------------------------------------------------------------------    
; select sky apertures and sky subtract
; ---------------------------------------------------------------------------    

;   iskyselect, skyfile=skyfile, skyapfile=skyapfile, norder_sky=1, /skyremember
    iall, paramfile, skyapfile=skyapfile

; ---------------------------------------------------------------------------    
; make the sensitivity curves - no standards were observed on Nov 12
; ---------------------------------------------------------------------------    

; generate sensitivity curves for each night (2.5")

    senstitle = '2001 '+['November 10','November 11','November 13']+'(2.5" Slit)'
    senslist = ['sens_2.5_01nov10.fits','sens_2.5_01nov11.fits','sens_2.5_01nov13.fits']
    iall, paramfile, stdfile=stdfile2_5, senstitle=senstitle, sensinfo=sensinfo, $
      slit_width=2.5, grey=2, doplot=doplot, sensname=senslist

;   readcol, stdfile2_5[i], stdlist, format='A', /silent, comment='#'
;   info = isensfunc(i1dnames(stdlist,aperture=21.67),grey=2,sensname=senslist[i],$
;     senstitle=senstitle[i],doplot=doplot,/write,slit_width=2.5,zptshift=0.0)

; generate sensitivity curves for each night (4.5")

    senstitle = '2001 '+['November 10','November 11','November 13']
    iall, paramfile, stdfile=stdfile4_5, senstitle=senstitle, sensinfo=sensinfo, $
      slit_width=4.5, grey=2, doplot=doplot

;   senslist = ['sens_4.5_01nov10.fits','sens_4.5_01nov11.fits','sens_4.5_01nov13.fits']
;   readcol, stdfile[i], stdlist, format='A', /silent, comment='#'
;   info = isensfunc(i1dnames(stdlist,aperture=21.67),grey=2,sensname=senslist[i],$
;     senstitle=senstitle[i],doplot=doplot,/write,slit_width=4.5,zptshift=0.0)

; make the combined sensitivity function 
    
    sensname = 'sens_2.5_01nov.fits' & senstitle='2001 November (2.5" Slit)'
    readcol, stdallfile2_5, stdlist, format='A', /silent, comment='#'
    ibatch, paramfile[0], stdlist=stdlist, sensname=sensname, /makesens, zptshift=0.0, $
      grey=2, sensinfo=sensinfo, senstitle=senstitle, slit_width=2.5, doplot=doplot, $
      wmapname='wmaplist_01nov.txt'

;   info = isensfunc(i1dnames(stdlist,aperture=21.67),grey=2,sensname=sensname,$
;     senstitle=senstitle,doplot=doplot,/write,slit_width=2.5,zptshift=0.0)

    sensname = 'sens_4.5_01nov.fits' & senstitle='2001 November (4.5" Slit)'
    readcol, stdallfile4_5, stdlist, format='A', /silent, comment='#'
    ibatch, paramfile[0], stdlist=stdlist, sensname=sensname, /makesens, zptshift=zptshift, $
      grey=2, sensinfo=sensinfo, senstitle=senstitle, slit_width=4.5, doplot=doplot, $
      wmapname='wmaplist_01nov.txt'

;   info = isensfunc(i1dnames(stdlist,aperture=21.67),grey=2,sensname=sensname,$
;     senstitle=senstitle,doplot=doplot,/write,slit_width=4.5,zptshift=zptshift)

; ---------------------------------------------------------------------------    
; generate the telluric spectra
; ---------------------------------------------------------------------------    

    tellindx = [0,1,3]
    for ii = 0L, n_elements(tellindx)-1L do begin

       i = tellindx[ii]

       readcol, tellfile[i], stdlist, format='A', /silent, comment='#'
       hotlist = i1dnames(stdlist,aperture=21.67)

       iconstruct_telluric, hotlist, tellmethod=2L, contmethod=1L, method=3L, $
         npoly=2, psname='qaplot_'+repstr(tellfits[i],'.fits','.ps'), tellfits=tellfits[i], $
         doplot=doplot, debug=0, /write

    endfor

; ---------------------------------------------------------------------------    
; flux-calibrate
; ---------------------------------------------------------------------------    
    
    sensname = 'sens_4.5_01nov.fits'
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

    label = [ ['10 November','Clear'],['11 November','Clear'],['13 November','Early clouds'],['2.5" Slit',''] ]
    senslist = ['sens_4.5_01nov10.fits','sens_4.5_01nov11.fits','sens_4.5_01nov13.fits','sens_2.5_01nov.fits']
    meansens = 'sens_4.5_01nov.fits'
    isenscompare, senslist, meansens=meansens, label=label, psname='qaplot_sens_compare_01nov.ps', /postscript

    loglist = ['qalog_sens_4.5_01nov10.log','qalog_sens_4.5_01nov11.log',$
      'qalog_sens_4.5_01nov13.log','qalog_sens_2.5_01nov.log','qalog_sens_4.5_01nov.log']
    label = [ ['10 November','Clear'],['11 November','Clear'],['13 November','Early clouds'],$
      ['2.5" Slit',''],['4.5" Slit',''] ]
    iqaphoto, loglist, label=label, title='2001 November', psname='qaplot_photo_01nov.ps', /postscript
    
; ---------------------------------------------------------------------------    
; generate the web page
; ---------------------------------------------------------------------------    

    if keyword_set(web) then begin
       readcol, 'weblist_01nov.txt', weblist, format='A', comment='#', /silent
       ispec_webpage, '01nov', weblist=weblist, html_path=atlas_path(/dataweb), $
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

; fix individual frames for cosmic rays

;   ibatch, paramfile[1], crlist='ra.3514.fits', /crclean, objlim=1.3, iaxis=1 ; N7137
;   ibatch, paramfile[1], crlist='ra.3524.fits', /crclean, iaxis=1 ; N337
;   ibatch, paramfile[1], crlist='ra.3520_2.fits', /crclean, objlim=1.3, sigfrac=1.0 ; N925W DRIFT55
;   ibatch, paramfile[1], crlist='ra.3532.fits', /crclean, objlim=1.0 ; N1377
    
;   ibatch, paramfile[2], crlist='ra.3721_2.fits', /crclean, objlim=1.5, iaxis=1 ; N24
;   ibatch, paramfile[2], crlist='ra.3758.fits', /crclean, objlim=1.8 ; M81 DWB
;   ibatch, paramfile[2], crlist='ra.3762_2.fits', /crclean, objlim=1.5 ; MRK33

; replace data 

;   rcube = rd2dspec('ra.3516.fits') ; not a cosmic ray jm02nov25uofa
;   ccube = rd2dspec('cra.3516.fits')
;   ccube.image[850:870,61:65] = rcube.image[850:870,61:65] & ccube.mask[850:870,61:65] = fix(0)
;   wrt2dspec, ccube.fname, ccube.image, ccube.sigmamap, ccube.mask, *ccube.header
;   icleanup, rcube & icleanup, ccube

;   ibatch, paramfile[1], crlist='ra.3534.fits', /crclean, iaxis=1

;   imcube = rd2dspec('ra.3534.fits') ; cosmic ray!
;   mask = imcube.mask*fix(0) & mask[845:885,36:42] = fix(32)
;   imcube.image = djs_maskinterp(imcube.image,mask,iaxis=1L)
;   imcube.mask = imcube.mask + mask
;   sxaddhist, 'Pixels [845:885,36:42] repaired '+im_today()+'.', *imcube.header
;   wrt2dspec, imcube.fname, imcube.image, imcube.sigmamap, imcube.mask, *imcube.header
;   icleanup, imcube

;   rcube = rd2dspec('ra.3534.fits') ; not a cosmic ray jm02nov25uofa
;   ccube = rd2dspec('cra.3534.fits')
;   ccube.image[951:955,62:66] = rcube.image[951:955,62:66] & ccube.mask[951:955,62:66] = fix(0)
;   wrt2dspec, ccube.fname, ccube.image, ccube.sigmamap, ccube.mask, *ccube.header
;   icleanup, rcube & icleanup, ccube

;   rcube = rd2dspec('ra.3550.fits') ; not a cosmic ray jm02nov25uofa
;   ccube = rd2dspec('cra.3550.fits')
;   ccube.image[702:717,85:88] = rcube.image[702:717,85:88] & ccube.mask[702:717,85:88] = fix(0)
;   wrt2dspec, ccube.fname, ccube.image, ccube.sigmamap, ccube.mask, *ccube.header
;   icleanup, rcube & icleanup, ccube

;   imcube = rd2dspec('ra.3744_2.fits') ; cosmic ray! ; jm02nov26uofa
;   mask = imcube.mask*fix(0) & mask[504:540,88:91] = fix(32)
;   imcube.image = djs_maskinterp(imcube.image,mask,iaxis=1L)
;   imcube.mask = imcube.mask + mask
;   sxaddhist, 'Pixels [504:540,88:91] repaired '+im_today()+'.', *imcube.header
;   wrt2dspec, imcube.fname, imcube.image, imcube.sigmamap, imcube.mask, *imcube.header
;   icleanup, imcube

;   ibatch, paramfile[2], crlist='ra.3744_2.fits', /crclean, objlim=1.0, iaxis=1
    
return
end
