pro mar04_script, doplot=doplot
; jm04mar12uofa

; ###########################################################################    
; THIS SCRIPT IS NOT COMPATIBLE WITH ISPECV2.0
; ###########################################################################    

stop    
    
    paramfile = ['ibatch_04mar16.txt']
    procfile = ['objlist_04mar16.txt']
    skyapfile = ['skyaplist_04mar16.txt']
    crsplitfile = ['crsplits_04mar16.txt']
    crfile = ['crlist_04mar16.txt']
;   tracefile = ['tracelist_04mar16.txt']
    stdfile = ['stdlist_4.5_04mar16.txt']
    stdfile2_5 = ['stdlist_2.5_04mar.txt']
    stdallfile = 'stdlist_4.5_04mar.txt'
    calibfile = ['caliblist_04mar16.txt']
    tellfile = ['telllist_04mar16.txt']

; ---------------------------------------------------------------------------
; pre-processing
; ---------------------------------------------------------------------------
    
;   iheader_check, root='a.'
;   ispec_makelists, ['a.04??'], ['04mar16'], /all, /overwrite, /gzip
    
    headfile = 'header_keywords_04mar.dat'
;   imake_keywordfile, ['SCANLEN','PA','APERTURE'], comments=['scan length [arcsec]',$
;     'slit position angle [degrees]','slit aperture [arcsec]'], values=[0.0,90.0,2.5], $
;     root='a.04??', headfile=headfile
    iheader_keywords, headfile, /update, silent=silent

; 04mar16
    
;   biases = file_search('a.0490[0-2]?.fits')
;   icalibavg, biases, outname='a.0401.fits', sigrej=3.0, /wfits
;   
;   domeflats = file_search('a.0490[5-7]?.fits')
;   icalibavg, domeflats, outname='a.0402.fits', sigrej=3.0, /wfits
;
;   skyflats = file_search('a.040[3-6].fits')
;   calib_combine, skyflats, 'skyflat_04mar16.fits', /sky, nhigh=1, nlow=1, $
;     overscan=[1208,3,1218,142]
    
; ---------------------------------------------------------------------------    
; carry out all the initial reductions
; ---------------------------------------------------------------------------    

    iall, paramfile, procfile=procfile, tracefile=tracefile, skyapfile=skyapfile, $
      crsplitfile=crsplitfile, crfile=crfile, doplot=doplot

; ---------------------------------------------------------------------------    
; make the sensitivity curves
; ---------------------------------------------------------------------------    

; generate sensitivity curves for each night (4.5")

    senstitle = '2004 '+['March 16']+' (4.5" Slit)'
    iall, paramfile, stdfile=stdfile, senstitle=senstitle, sensinfo=sensinfo, grey=2, $
      slit_width=4.5, doplot=doplot, /telluric

;   senslist = ['sens_4.5_04mar16.fits']
;   readcol, stdfile[i], stdlist, format='A', /silent, comment='#'
;   info = isensfunc(i1dnames(stdlist,aperture=22.5),grey=2,sensname=senslist[i],$
;     senstitle=senstitle[i],/doplot,/write,slit_width=4.5,/telluric)

; make the combined sensitivity function
    
    sensname = 'sens_2.5_04mar.fits' & senstitle = '2004 March (2.5" Slit)'
    readcol, stdfile2_5, stdlist, format='A', /silent, comment='#'
    ibatch, paramfile[0], stdlist=stdlist, sensname=sensname, /makesens, grey=2, $
      sensinfo=sensinfo, senstitle=senstitle, slit_width=2.5, doplot=doplot, $
      wmapname='wmaplist_04mar.txt', /telluric

;   info = isensfunc(i1dnames(stdlist,aperture=21),grey=2,sensname=sensname,$
;     senstitle=senstitle,doplot=doplot,/write,slit_width=2.5,zptshift=0.0,/telluric)

    sensname = 'sens_4.5_04mar.fits' & senstitle='2004 March (4.5" Slit)'
    readcol, stdallfile, stdlist, format='A', /silent, comment='#'
    ibatch, paramfile[0], stdlist=stdlist, sensname=sensname, /makesens, grey=2, $
      sensinfo=sensinfo, senstitle=senstitle, slit_width=4.5, doplot=doplot, $
      wmapname='wmaplist_04mar.txt', /telluric

;   info = isensfunc(i1dnames(stdlist,aperture=22.5),grey=2,sensname=sensname,$
;     senstitle=senstitle,/doplot,/write,slit_width=4.5,/telluric)

; ---------------------------------------------------------------------------    
; calibrate the data and append the nightly telluric absorption spectrum
; ---------------------------------------------------------------------------    
    
    sensname = 'sens_4.5_04mar.fits'
    iall, paramfile, calibfile=calibfile, tellfile=tellfile, sensname=sensname

; ---------------------------------------------------------------------------    
; miscellaneous tasks
; ---------------------------------------------------------------------------    

; compute the seeing

    iseeing, 'seeing_04mar.txt', title='2004 March Seeing', psname='seeing_04mar.ps', $
      /postscript, doplot=doplot

; compare the sensitivity curves and compute the quality assurance on
; the photometry
;
;   label = reform(['21 March','2.5" Slit',''],3,1)
;   senslist = 'sens_2.5_04mar14.fits'
;   meansens = 'sens_4.5_04mar14.fits'
;   isenscompare, senslist, meansens=meansens, title='1998 March 21', $
;     label=label, psname='qaplot_sens_compare_98mar.ps', /postscript
;   
;   loglist = ['qalog_sens_2.5_04mar14.log','qalog_sens_4.5_04mar14.log']
;   label = [ ['21 March','2.5" Slit','Clear'],['21 March','4.5" Slit','Clear'] ]
;   iqaphoto, loglist, label=label, title='1998 March 21', psname='qaplot_photo_04mar14.ps', /postscript

; ---------------------------------------------------------------------------    
; generate the web page
; ---------------------------------------------------------------------------    

;   readcol, 'weblist_04mar.txt', weblist, format='A', comment='#', /silent
;   ispec_webpage, '04mar', rootname='fw', weblist=weblist, red_path=red_path, /cleanpng, $
;     html_path=atlas_path(/dataweb), html_only=html_only

return
end
