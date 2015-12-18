pro feb02_script, doplot=doplot, zptshift=zptshift, web=web

    dateroot = '02feb'
    dates = dateroot+['07','08','09','10']

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
    
;   crsplitfile = ['crsplits_02feb07.txt','crsplits_02feb08.txt','','crsplits_02feb10.txt']

; ---------------------------------------------------------------------------
; pre-processing
; ---------------------------------------------------------------------------
    
;   iheader_check, root='a.'
;   ispec_makelists, ['a.38','a.39','a.40','a.41'], ['02feb07','02feb08','02feb09','02feb10'], $
;     /all, /overwrite, /gzip
    
;   skyflats = file_search('a.'+['380[2-4]','390[3-5]','400[3-5]','410[3-5]']+'.fits')
;   skyflatfile = 'skyflat_02feb.fits'
;   arm_flatcombine, skyflats, skyflatfile, [1206,3,1218,122], /display, psfile='qaplot_skyflat_02feb.ps'
    
;   imake_keywordfile, ['SCANLEN','PA','APERTURE'], comments=['scan length [arcsec]',$
;     'slit position angle [degrees]','slit aperture [arcsec]'], values=[0.0,90.0,2.5], $
;     root='a.01', headfile=headfile

; correct some headers

;   list = ['a.3850.fits','a.3852.fits','a.3853.fits','a.3854.fits']
;   object = ['N4450 nuc','N4450 drift20','N4450 drift55','N4450 drift55']
;   nlist = n_elements(list)
;   for i = 0L, nlist-1L do begin
;      h = headfits(list[i])
;      sxaddpar, h, 'object', object[i]
;      modfits, list[i], 0, h
;   endfor

; update the headers    
    
;   headfile = 'header_keywords_02feb.dat'
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

; NGC3184 East drift 
    icrcombine, 'ra.'+['4124','4125']+'.fits', /find_pixshift, $
      user_pixshift=[0.0,2.0], /wfits;, /debug
    ibatch, paramfile[3], crlist='ra.4124_2.fits', /crclean

; NGC3184 West drift 
    icrcombine, 'ra.'+['4126','4127']+'.fits', /find_pixshift, $
      user_pixshift=[0.0,0.0], /wfits;, /debug
    ibatch, paramfile[3], crlist='ra.4126_2.fits', /crclean

; NGC3351 West drift 
    icrcombine, 'ra.'+['4132','4133']+'.fits', /find_pixshift, $
      user_pixshift=[0.0,0.0], /wfits;, /debug
    ibatch, paramfile[3], crlist='ra.4132_2.fits', /crclean

; clean up residual cosmic rays; 3930 isn't perfect but that's okay 

    ibatch, paramfile[0], crlist='ra.3856.fits', /crclean, objlim=1.0, iaxis=0

; repair LA_COSMIC pixels

    repair_la_cosmic_crpix, ['cra.3944.fits','cra.3945.fits'], $
      ['repair.3944.dat','repair.3945.dat'], /repair
    repair_la_cosmic_crpix, 'cra.3925_2.fits', 'repair.3925_2.dat', /repair

; ---------------------------------------------------------------------------    
; SINGS: stitch; select sky apertures and sky-subtract
; ---------------------------------------------------------------------------    

    sings2d_stitch, 'cra.'+['4124_2','4126_2']+'.fits', outname='wcra.4124_2_stitch.fits', $
      object='NGC3184 drift55', user_pixshift=[99.0], wmapname='wmap_ra.4134.idlsave', $
      refindx=0L, minwave=3635.0, dwave=2.75, scalefactor=[1.0,1.04], debug=1, wfits=1

    sings2d_stitch, 'cra.'+['4130_2','4132_2']+'.fits', outname='wcra.4130_2_stitch.fits', $
      object='NGC3351 drift55', user_pixshift=[104.0], wmapname='wmap_ra.4134.idlsave', $
      refindx=0L, minwave=3635.0, dwave=2.75, scalefactor=[1.02,1.0], debug=1, wfits=1

;   readcol, 'skylist_stitch.txt', skylist, format='A', comment='#', /silent
;   iskyselect, skylist=skylist, skyapfile='skyaplist_stitch.txt'
    ibatch, paramfile[3], skyapfile='skyaplist_stitch.txt', skymethod=3L, /skysub
    
; ---------------------------------------------------------------------------    
; select sky apertures and sky subtract
; ---------------------------------------------------------------------------    

;   iskyselect, skyfile=skyfile, skyapfile=skyapfile, norder_sky=1, /skyremember

    iall, paramfile, skyapfile=skyapfile

; ---------------------------------------------------------------------------    
; make the sensitivity curves
; ---------------------------------------------------------------------------    

; generate sensitivity curves for each night (2.5")

    senstitle = '2002 '+['February 07','February 08','February 09','February 10']+'( 2.5" Slit)'
    senslist = ['sens_2.5_02feb07.fits','sens_2.5_02feb08.fits','sens_2.5_02feb09.fits','sens_2.5_02feb10.fits']
    iall, paramfile, stdfile=stdfile2_5, senstitle=senstitle, sensinfo=sensinfo, $
      slit_width=2.5, grey=2, doplot=doplot, sensname=senslist

;   readcol, stdfile2_5[i], stdlist, format='A', /silent, comment='#'
;   info = isensfunc(i1dnames(stdlist,aperture=20),grey=2,sensname=senslist[i],$
;     senstitle=senstitle[i],doplot=doplot,/write,slit_width=2.5,zptshift=0.0)

; generate sensitivity curves for each night (4.5")

    senstitle = '2002 '+['February 07','February 08','February 09','February 10']
    iall, paramfile, stdfile=stdfile4_5, senstitle=senstitle, sensinfo=sensinfo, $
      slit_width=4.5, grey=2, doplot=doplot

;   senslist = ['sens_4.5_02feb07.fits','sens_4.5_02feb08.fits','sens_4.5_02feb09.fits','sens_4.5_02feb10.fits']
;   readcol, stdfile[i], stdlist, format='A', /silent, comment='#'
;   info = isensfunc(i1dnames(stdlist,aperture=20),grey=2,sensname=senslist[i],$
;     senstitle=senstitle[i],doplot=doplot,/write,slit_width=4.5,zptshift=0.0)

; make the combined sensitivity functions
    
    sensname = 'sens_2.5_02feb.fits' & senstitle='2002 February (2.5" Slit)'
    readcol, stdallfile2_5, stdlist, format='A', /silent, comment='#'
    ibatch, paramfile[0], stdlist=stdlist, sensname=sensname, /makesens, zptshift=0.0, $
      grey=2, sensinfo=sensinfo, senstitle=senstitle, slit_width=2.5, doplot=doplot, $
      wmapname='wmaplist_02feb.txt'

;   info = isensfunc(i1dnames(stdlist,aperture=20),grey=2,sensname=sensname,$
;     senstitle=senstitle,doplot=doplot,/write,slit_width=2.5,zptshift=0.0)

    sensname = 'sens_4.5_02feb.fits' & senstitle='2002 February (4.5" Slit)'
    readcol, stdallfile4_5, stdlist, format='A', /silent, comment='#'
    ibatch, paramfile[0], stdlist=stdlist, sensname=sensname, /makesens, zptshift=zptshift, $
      grey=2, sensinfo=sensinfo, senstitle=senstitle, slit_width=4.5, doplot=doplot, $
      wmapname='wmaplist_02feb.txt'

;   info = isensfunc(i1dnames(stdlist,aperture=20),grey=2,sensname=sensname,$
;     senstitle=senstitle,doplot=doplot,/write,slit_width=4.5,zptshift=zptshift)

; ---------------------------------------------------------------------------    
; generate the telluric spectra
; ---------------------------------------------------------------------------    

    for i = 0L, n_elements(paramfile)-1L do begin

       readcol, tellfile[i], stdlist, format='A', /silent, comment='#'
       hotlist = i1dnames(stdlist,aperture=20)

       iconstruct_telluric, hotlist, tellmethod=2L, contmethod=1L, method=3L, $
         npoly=2, psname='qaplot_'+repstr(tellfits[i],'.fits','.ps'), tellfits=tellfits[i], $
         doplot=doplot, debug=0, /write

    endfor
    
; ---------------------------------------------------------------------------    
; flux-calibrate
; ---------------------------------------------------------------------------    
    
    sensname = 'sens_4.5_02feb.fits'
    iall, paramfile, fluxcalibfile=calibfile, calibfile=calibfile, $
      sensname=sensname;, tellfits=tellfits

; the stitched spectra have been wavelength calibrated already, so
; just flux calibrate    
    
    readcol, 'skylist_stitch.txt', caliblist, format='A', comment='#', /silent
    ibatch, paramfile[3], caliblist='s'+caliblist, $ ; tellfits=tellfits[3], $
      sensname=sensname, /fluxcalibrate

; ---------------------------------------------------------------------------    
; miscellaneous tasks
; ---------------------------------------------------------------------------    

; compare the sensitivity curves

    label = [ ['07 February','Clear'],['08 February','Clear'],['09 February','Clear, windy'],$
      ['10 February','Clear, windy'],['2.5" Slit',''] ]
    senslist = ['sens_4.5_02feb07.fits','sens_4.5_02feb08.fits','sens_4.5_02feb09.fits',$
      'sens_4.5_02feb10.fits','sens_2.5_02feb.fits']
    meansens = 'sens_4.5_02feb.fits'
    isenscompare, senslist, meansens=meansens, label=label, title='2002 February', $
      psname='qaplot_sens_compare_02feb.ps', /postscript

    loglist = ['qalog_sens_4.5_02feb07.log','qalog_sens_4.5_02feb08.log',$
      'qalog_sens_4.5_02feb09.log','qalog_sens_4.5_02feb10.log','qalog_sens_2.5_02feb.log',$
      'qalog_sens_4.5_02feb.log']
    label = [ ['07 February','Clear'],['08 February','Clear'],['09 February','Clear, windy'],$
      ['10 February','Clear, windy'],['2.5" Slit',''],['4.5" Slit',''] ]
    iqaphoto, loglist, label=label, title='2002 February', psname='qaplot_photo_02feb.ps', /postscript

; ---------------------------------------------------------------------------    
; generate the web page
; ---------------------------------------------------------------------------    

    if keyword_set(web) then begin
       readcol, 'weblist_02feb.txt', weblist, format='A', comment='#', /silent
       ispec_webpage, '02feb', weblist=weblist, html_path=atlas_path(/dataweb), $
         html_only=html_only, /cleanpng
    endif

; ---------------------------------------------------------------------------    
; old code
; ---------------------------------------------------------------------------    

; fix individual frames for cosmic rays

;   ibatch, paramfile[0], crlist='ra.3856.fits', /crclean, objlim=1.5 ; N5194
;   ibatch, paramfile[1], crlist='ra.3933.fits', /crclean, objlim=1.0, sigclip=4.0 ; N3265
;   ibatch, paramfile[2], crlist='ra.4023_2.fits', /crclean, objlim=1.5 ; N3351
;   ibatch, paramfile[2], crlist='ra.4036.fits', /crclean, objlim=1.5 ; N3239
;
;   imcube = rd2dspec('ra.4130_2.fits') ; bad row induced by a bad cosmic ray
;   mask = imcube.mask*fix(0) & mask[1139:1152,63:64] = fix(32)
;   imcube.image = djs_maskinterp(imcube.image,mask,iaxis=1L)
;   imcube.mask = imcube.mask + mask
;   sxaddhist, 'Pixels [1139:1152,63:64] repaired '+im_today()+'.', *imcube.header
;   wrt2dspec, imcube.fname, imcube.image, imcube.sigmamap, imcube.mask, *imcube.header
;   icleanup, imcube
;   ibatch, paramfile[3], crlist='ra.4130_2.fits', /crclean
    
return
end
