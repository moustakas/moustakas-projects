pro jun02_script, doplot=doplot, zptshift=zptshift
    
    if n_elements(zptshift) eq 0L then zptshift = 0.0 ; sensitivity function zero point shift (mag)
    
    paramfile = ['ibatch_02jun25.txt','ibatch_02jun26.txt']
    procfile = ['objlist_02jun25.txt','objlist_02jun26.txt']
    skyfile = ['skylist_02jun25.txt','skylist_02jun26.txt']
    skyapfile = ['skyaplist_02jun25.txt','skyaplist_02jun26.txt']
    crsplitfile = ['','crsplits_02jun26.txt']
    crfile = ['crlist_02jun25.txt','crlist_02jun26.txt']
    stdfile = ['stdlist_4.5_02jun25.txt','stdlist_4.5_02jun26.txt']
    stdallfile = 'stdlist_4.5_02jun.txt'
    calibfile = ['caliblist_02jun25.txt','caliblist_02jun26.txt']
    tellfile = ['telllist_02jun25.txt','telllist_02jun26.txt']

; ---------------------------------------------------------------------------
; pre-processing
; ---------------------------------------------------------------------------
    
;   iheader_check, root='a.'
;   ispec_makelists, ['a.48','a.49'], ['02jun25','02jun26'], /all, /overwrite, /gzip
    
;   imake_keywordfile, ['SCANLEN','PA','APERTURE'], comments=['scan length [arcsec]',$
;     'slit position angle [degrees]','slit aperture [arcsec]'], values=[0.0,90.0,2.5], $
;     root='a.01', headfile=headfile

    headfile = 'header_keywords_02jun.dat'
    iheader_keywords, headfile, /update, silent=silent

; ---------------------------------------------------------------------------
; correct headers
; ---------------------------------------------------------------------------
    
    h = headfits('a.4914.fits')
    sxaddpar, h, 'OBJECT', 'HD167946'
    modfits, 'a.4914.fits', 0, h
    
; ---------------------------------------------------------------------------
; carry out all the initial reductions
; ---------------------------------------------------------------------------
    
    iall, paramfile, procfile=procfile, tracefile=tracefile, doplot=doplot

;   iall, paramfile, procfile=procfile, tracefile=tracefile, skyfile=skyfile, $
;     skyapfile=skyapfile, crsplitfile=crsplitfile, crfile=crfile, doplot=doplot
    
; ---------------------------------------------------------------------------    
; select sky apertures    
; ---------------------------------------------------------------------------    
   
;   iskyselect, skyfile=skyfile, skyapfile=skyapfile, norder_sky=1, /skyremember

; ---------------------------------------------------------------------------    
; sky subtract, combine cr-splits, and cosmic-ray reject
; ---------------------------------------------------------------------------    

    iall, paramfile, skyapfile=skyapfile, crsplitfile=crsplitfile, $
      crfile=crfile, doplot=doplot

; ---------------------------------------------------------------------------    
; make the sensitivity curves
; ---------------------------------------------------------------------------    

; generate sensitivity curves for each night (4.5")

    senstitle = '2002 '+['June 25','June 26']
    iall, paramfile, stdfile=stdfile, senstitle=senstitle, sensinfo=sensinfo, $
      slit_width=4.5, grey=2, doplot=doplot;, /telluric

;   senslist = ['sens_4.5_02jun25.fits','sens_4.5_02jun26.fits']
;   readcol, stdfile[i], stdlist, format='A', /silent, comment='#'
;   info = isensfunc(i1dnames(stdlist,aperture=22.17),grey=2,sensname=senslist[i],$
;     senstitle=senstitle[i],doplot=doplot,/write,slit_width=4.5,zptshift=0.0);,/telluric)

; make the combined sensitivity function

    sensname = 'sens_4.5_02jun.fits' & senstitle='2002 June (4.5" Slit)'
    readcol, stdallfile, stdlist, format='A', /silent, comment='#'
    ibatch, paramfile[0], stdlist=stdlist, sensname=sensname, /makesens, zptshift=zptshift, $
      grey=2, sensinfo=sensinfo, senstitle=senstitle, slit_width=4.5, doplot=doplot, $
      wmapname='wmaplist_02jun.txt';, /telluric

;   info = isensfunc(i1dnames(stdlist,aperture=22.17),grey=2,sensname=sensname,$
;     senstitle=senstitle,doplot=doplot,/write,slit_width=4.5,zptshift=zptshift);,/telluric)

; ---------------------------------------------------------------------------    
; flux-calibrate
; ---------------------------------------------------------------------------    

    sensname = 'sens_4.5_02jun.fits'
    iall, paramfile, calibfile=calibfile, sensname=sensname
;   iall, paramfile, calibfile=calibfile, tellfile=tellfile, sensname=sensname

; ---------------------------------------------------------------------------    
; miscellaneous tasks
; ---------------------------------------------------------------------------    

; compute the seeing

    iseeing, 'seeing_02jun.txt', title='2002 June Seeing', psname='seeing_02jun.ps', /postscript, doplot=doplot

; compare the sensitivity curves

    label = [ ['25 June','Cloudy'],['26 June','Cloudy'] ]
    senslist = ['sens_4.5_02jun25.fits','sens_4.5_02jun26.fits']
    meansens = 'sens_4.5_02jun.fits'
    isenscompare, senslist, meansens=meansens, title='2002 June', label=label, $
      psname='qaplot_sens_compare_02jun.ps', /postscript
    
    loglist = ['qalog_sens_4.5_02jun25.log','qalog_sens_4.5_02jun26.log','qalog_sens_4.5_02jun.log']
    label = [ ['25 June','Cloudy'],['26 June','Cloudy'],['4.5" slit',''] ]
    iqaphoto, loglist, label=label, title='2002 June', psname='qaplot_photo_02jun.ps', /postscript
    
; ---------------------------------------------------------------------------    
; generate the web page
; ---------------------------------------------------------------------------    

    readcol, 'weblist_02jun.txt', weblist, format='A', comment='#', /silent
    ispec_webpage, '02jun', rootname='fw', weblist=weblist, red_path=red_path, /cleanpng, $
      html_path=atlas_path(/dataweb), html_only=html_only

; ---------------------------------------------------------------------------    
; old code
; ---------------------------------------------------------------------------    
    
; repair specific images for cosmic rays

;   ibatch, paramfile[1], crlist='ra.4915.fits', /crclean, sigfrac=1.0, objlim=2.5 ; HD166384

return
end
