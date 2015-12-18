pro may99_script, doplot=doplot, zptshift=zptshift, web=web
; jm02jan3uofa

    dateroot = '99may'
    dates = dateroot+['08','09','10','11']

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
;   ispec_makelists, ['a.20','a.21','a.22','a.23'], ['99may08','99may09','99may10','99may11'], $
;     /all, /overwrite, gzip=0
    
;   skyflats = file_search('a.'+['200[4-6]','220[3-4]','230[4-5]']+'.fits')
;   skyflatfile = 'skyflat_99may.fits'
;   arm_flatcombine, skyflats, skyflatfile, [1208,5,1217,130], /display, psfile='qaplot_skyflat_99may.ps'
    
;   imake_keywordfile, ['SCANLEN','POSANGLE','APERTURE'], comments=['scan length [arcsec]',$
;     'slit position angle [degrees]','slit aperture [arcsec]'], values=[0.0,90.0,2.5], $
;     root='a.01', headfile=headfile

; correct some headers

;   h = headfits('a.2322.fits')
;   sxaddpar, h, 'OBJECT', 'DD0 165 drift'
;   modfits, 'a.2322.fits', 0, h
;   
;   h = headfits('a.2323.fits')
;   sxaddpar, h, 'OBJECT', 'DD0 190 drift'
;   modfits, 'a.2323.fits', 0, h

; update the headers    

;   headfile = 'header_keywords_99may.dat'
;   iheader_keywords, headfile, /update, silent=silent

; ---------------------------------------------------------------------------
; carry out all the initial reductions
; ---------------------------------------------------------------------------

    iall, paramfile, procfile=procfile, doplot=doplot, tracefile=tracefile

; ---------------------------------------------------------------------------    
; combine cr-splits and cosmic-ray reject
; ---------------------------------------------------------------------------    

; no pixel shifts are needed since there is only one crsplit file and
; for that object (IC0860) the pixel shift is zero    
    
;   ibatch, paramfile[2], crsplitfile=crsplitfile[2], /find_pixshift, /debug, /crsplits

    iall, paramfile, crsplitfile=crsplitfile, crfile=crfile;, /find_pixshift ; <-- NOTE!

; clean up residual cosmic rays

    ibatch, paramfile[1], crlist='ra.2111.fits', /crclean, objlim=1.0, sigclip=8.0, iaxis=0
    ibatch, paramfile[3], crlist='ra.2311.fits', /crclean, objlim=0.4, sigclip=13.0, iaxis=0
    ibatch, paramfile[3], crlist='ra.2313.fits', /crclean, objlim=0.5, sigclip=13.0, iaxis=0

; repair LA_COSMIC pixels

    repair_la_cosmic_crpix, ['cra.2109.fits','cra.2311.fits'], $
      ['repair.2109.dat','repair.2311.dat'], /repair

; ---------------------------------------------------------------------------    
; ATLAS: stitch
; ---------------------------------------------------------------------------    

; NGC3344    
    
    sings2d_stitch, 'cra.21'+['09','08']+'.fits', outname='wcra.2108_stitch.fits', $
      object='NGC3344 drift', scalefactor=[1.0,1.03], $
      user_pixshift=[107.0], wmapname='wmap_ra.2107.idlsave', $
      refindx=1L, minwave=3680.0, dwave=2.75, debug=1, wfits=1
    
;   spec1list = ['fscra.2108.fits','fwscra.2108.fits']
;   spec2list = ['fscra.2109.fits','fwscra.2109.fits']
;   stitchlist = ['fscra.2108-09.fits','fwscra.2108-09.fits']
;   objectlist = ['N3344 drift','N3344 drift']
;
;   atlas2d_stitch, spec1list, spec2list, stitchlist, objectlist, $
;     user_pixshift=[-19L,-19L], debug=debug, /wfits
    
; stitch together the NGC5194 East-West scan
 
    sings2d_stitch, 'cra.22'+['20','19']+'.fits', outname='wcra.2219_stitch.fits', $
      object='NGC5194 drift', scalefactor=[1.06,1.0], $
      user_pixshift=[109.0], wmapname='wmap_ra.2231.idlsave', $
      refindx=1L, minwave=3680.0, dwave=2.75, debug=1, wfits=1

;   spec1list = ['fscra.2219.fits','fwscra.2219.fits']
;   spec2list = ['fscra.2220.fits','fwscra.2220.fits']
;   stitchlist = ['fscra.2219-20.fits','fwscra.2219-20.fits']
;   objectlist = ['N5194 drift','N5194 drift']
;
;   atlas2d_stitch, spec1list, spec2list, stitchlist, objectlist, $
;     user_pixshift=[-18L,-18L], debug=debug, /wfits
    
;   stitch_ngc_5194, /wfits
;   stitch_ngc_3344, /wfits

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
    
; generate sensitivity curves for all the nights (2.5")

    senslist = ['sens_2.5_99may08.fits','sens_2.5_99may09.fits','sens_2.5_99may10.fits','sens_2.5_99may11.fits']
    senstitle = '1999 '+['May 08','May 09','May 10','May 11']+' (2.5" Slit)'
    iall, paramfile, stdfile=stdfile2_5, senstitle=senstitle, sensinfo=sensinfo, $
      grey=2, slit_width=2.5, doplot=doplot, sensname=senslist
    
;   readcol, stdfile[i+4], stdlist, format='A', /silent, comment='#'
;   info = isensfunc(i1dnames(stdlist,aperture=21),grey=2,sensname=senslist[i],$
;     senstitle=senstitle[i],doplot=doplot,/write,slit_width=2.5,zptshift=0.0)

; generate sensitivity curves for all the nights (4.5")

    senstitle = '1999 '+['May 08','May 09','May 10','May 11']
    iall, paramfile, stdfile=stdfile4_5, senstitle=senstitle, sensinfo=sensinfo, $
      grey=2, slit_width=4.5, doplot=doplot
    
;   senslist = ['sens_4.5_99may08.fits','sens_4.5_99may09.fits','sens_4.5_99may10.fits','sens_4.5_99may11.fits']
;   readcol, stdfile[i+4], stdlist, format='A', /silent, comment='#'
;   info = isensfunc(i1dnames(stdlist,aperture=21),grey=2,sensname=senslist[i],$
;     senstitle=senstitle[i],doplot=doplot,/write,slit_width=4.5,zptshift=0.0)

; generate a sensitivity curve for the drift-scanned star for May

    ibatch, paramfile[1], stdlist='scra.2124.fits', sensname='sens_drift_99may.fits', $
      sensinfo=sensinfo, senstitle='1999 May (Drift)', /makesens, slit_width=20.0, $
      doplot=doplot
    
; make the combined sensitivity function for the May run
    
    sensname = 'sens_2.5_99may.fits' & senstitle='1999 May (2.5" Slit)'
    readcol, stdallfile2_5, stdlist, format='A', /silent, comment='#'
    ibatch, paramfile[0], stdlist=stdlist, sensname=sensname, /makesens, zptshift=0.0, $
      grey=2, sensinfo=sensinfo, senstitle=senstitle, slit_width=2.5, doplot=doplot, $
      wmapname='wmaplist_99may.txt'

;   info = isensfunc(i1dnames(stdlist,aperture=21),grey=2,sensname=sensname,$
;     senstitle=senstitle,doplot=doplot,/write,slit_width=2.5,zptshift=0.0)

    sensname = 'sens_4.5_99may.fits' & senstitle='1999 May (4.5" Slit)'
    readcol, stdallfile4_5, stdlist, format='A', /silent, comment='#'
    ibatch, paramfile[0], stdlist=stdlist, sensname=sensname, /makesens, zptshift=zptshift, $
      grey=2, sensinfo=sensinfo, senstitle=senstitle, slit_width=4.5, doplot=doplot, $
      wmapname='wmaplist_99may.txt'

;   info = isensfunc(i1dnames(stdlist,aperture=21),grey=2,sensname=sensname,$
;     senstitle=senstitle,doplot=doplot,/write,slit_width=4.5,zptshift=zptshift)

; ---------------------------------------------------------------------------    
; generate the telluric spectra
; ---------------------------------------------------------------------------    

    for i = 0L, n_elements(paramfile)-1L do begin

       readcol, tellfile[i], stdlist, format='A', /silent, comment='#'
       hotlist = i1dnames(stdlist,aperture=21)

       iconstruct_telluric, hotlist, tellmethod=2L, contmethod=1L, method=3L, $
         npoly=2, psname='qaplot_'+repstr(tellfits[i],'.fits','.ps'), tellfits=tellfits[i], $
         doplot=doplot, debug=0, /write

    endfor
    
; ---------------------------------------------------------------------------    
; flux-calibrate
; ---------------------------------------------------------------------------    
    
    sensname = 'sens_4.5_99may.fits'
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

    label = [ ['08 May','Clear, Windy'],['09 May','Clear, Windy'],['10 May','Thin clouds'],$
      ['11 May','Clear'],['Drift',''],['2.5" Slit',''] ]
    senslist = ['sens_4.5_99may08.fits','sens_4.5_99may09.fits','sens_4.5_99may10.fits',$
      'sens_4.5_99may11.fits','sens_drift_99may.fits','sens_2.5_99may.fits']
    meansens = 'sens_4.5_99may.fits'
    isenscompare, senslist, meansens=meansens, title='1999 May', label=label, psname='qaplot_sens_compare_99may.ps', /postscript

    loglist = ['qalog_sens_4.5_99may08.log','qalog_sens_4.5_99may09.log','qalog_sens_4.5_99may10.log',$
      'qalog_sens_4.5_99may11.log','qalog_sens_drift_99may.log','qalog_sens_2.5_99may.log','qalog_sens_4.5_99may.log']
    label = [ ['08 May','Clear, Windy'],['09 May','Clear, Windy'],['10 May','Thin clouds'],$
      ['11 May','Clear'],['Drift',''],['2.5" Slit',''],['4.5" Slit',''] ]
    iqaphoto, loglist, label=label, title='1999 May', psname='qaplot_photo_99may.ps', /postscript
    
; ---------------------------------------------------------------------------    
; generate the web pages
; ---------------------------------------------------------------------------    

    if keyword_set(web) then begin
       readcol, 'weblist_99may.txt', weblist, format='A', comment='#', /silent
       ispec_webpage, '99may', weblist=weblist, html_path=atlas_path(/dataweb), $
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

; fix specific images for cosmic rays

;   ibatch, paramfile[0], crlist='ra.1617.fits', /crclean, objlim=1.3 ; A1367 pos 7
;   ibatch, paramfile[0], crlist='ra.1624.fits', /crclean, objlim=1.0, iaxis=1 ; N4656
;   ibatch, paramfile[0], crlist='ra.1627.fits', /crclean, objlim=1.5 ; N5430
;   ibatch, paramfile[0], crlist='ra.1629.fits', /crclean, objlim=1.3 ; N5953/4
;   ibatch, paramfile[5], crlist='ra.2111.fits', /crclean, objlim=1.0 ; U06456
;   ibatch, paramfile[5], crlist='ra.2113.fits', /crclean, objlim=1.5 ; N4144
;   ibatch, paramfile[6], crlist='ra.2211.fits', /crclean, objlim=1.2 ; N2500
;   ibatch, paramfile[6], crlist='ra.2215.fits', /crclean, objlim=1.5 ; U07690
;   ibatch, paramfile[6], crlist='ra.2225.fits', /crclean, objlim=1.5 ; IZw107
;   ibatch, paramfile[7], crlist='ra.2311.fits', /crclean, sigfrac=1.7 ; N3432
;   ibatch, paramfile[7], crlist='ra.2313.fits', /crclean, sigfrac=1.7 ; N3104
;   ibatch, paramfile[7], crlist='ra.2315.fits', /crclean, objlim=1.3 ; DDO161
;   ibatch, paramfile[7], crlist='ra.2323.fits', /crclean, objlim=1.3 ; DD0190

return
end
