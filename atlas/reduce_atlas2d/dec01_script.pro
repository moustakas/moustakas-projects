pro dec01_script, doplot=doplot, zptshift=zptshift, web=web

    zptshift = 0.0 ; sensitivity function zero point shift (mag)    

    dateroot = '01dec'
    dates = dateroot+['20','21','22','23','24']

    paramfile = 'ibatch_'+dates+'.txt'
    procfile = 'objlist_'+dates+'.txt'
    skyfile = 'skylist_'+dates+'.txt'
    skyapfile = 'skyaplist_'+dates+'.txt'
    tracefile = 'tracelist_'+dates+'.txt'
    crsplitfile = 'crsplits_'+dates+'.txt'
    crfile = 'crlist_'+dates+'.txt'
    stdallfile = 'stdlist_'+dateroot+'.txt'
    calibfile = 'caliblist_'+dates+'.txt'
    tellfile = 'stdlist_'+dates+'.txt'
    tellfits = 'telluric_'+dates+'.fits'
    
; ---------------------------------------------------------------------------
; pre-processing
; ---------------------------------------------------------------------------

;   iheader_check, root='n[1-5]'
;   ispec_makelists, ['n1','n2','n3','n5'], ['01dec20','01dec21','01dec22','01dec24'], /all, /overwrite, /gzip
    
;   imake_keywordfile, ['SCANLEN','POSANGLE','APERTURE'], comments=['scan length [arcsec]',$
;     'slit position angle [degrees]','slit aperture [arcsec]'], values=[0.0,90.0,3.0], $
;     root='n?', headfile=headfile

; Dec 20
;
;   biasfile = 'n1.zero.fits'
;   biases = file_search(['n100[0-2]?.fits','n10030.fits'])
;   arm_zerocombine, biases, biasfile, display=doplot, psfile='qaplot_bias_01dec20.ps'
;
;   domefile = 'n1.domeflat.fits'
;   domeflats = file_search(['n1003[2-9].fits','n100[4-5]?.fits','n10060.fits'])
;   arm_flatcombine, domeflats, domefile, [1225,50,1265,286], $
;     display=doplot, psfile='qaplot_domeflat_01dec20.ps'
;
;   skyflatfile = 'n1.skyflat.fits'
;   skyflats = file_search('n1006[3-8].fits')
;   arm_flatcombine, skyflats, skyflatfile, [1225,50,1265,286], $
;     display=doplot, psfile='qaplot_skyflat_01dec20.ps'
;
; Dec 21
;
;   biasfile = 'n2.zero.fits'
;   biases = file_search(['n200[0-2]?.fits','n20030.fits'])
;   arm_zerocombine, biases, biasfile, display=doplot, psfile='qaplot_bias_01dec21.ps'
;
;   domefile = 'n2.domeflat.fits'
;   domeflats = file_search(['n2003[1-9].fits','n200[4-5]?.fits','n20060.fits'])
;   arm_flatcombine, domeflats, domefile, [1225,50,1265,286], $
;     display=doplot, psfile='qaplot_domeflat_01dec21.ps'
;
;   skyflatfile = 'n2.skyflat.fits'
;   skyflats = file_search(['n2006[6-9].fits','n2007[0-1].fits'])
;   arm_flatcombine, skyflats, skyflatfile, [1225,50,1265,286], $
;     display=doplot, psfile='qaplot_skyflat_01dec21.ps'
;
; Dec 22
;
;   biasfile = 'n3.zero.fits'
;   biases = file_search(['n300[0-2]?.fits','n30030.fits'])
;   arm_zerocombine, biases, biasfile, display=doplot, psfile='qaplot_bias_01dec22.ps'
;
;   domefile = 'n3.domeflat.fits'
;   domeflats = file_search(['n3003[1-9].fits','n3004[0-5].fits',$
;     'n3004[7-9].fits','n3005?.fits','n30060.fits'])
;   arm_flatcombine, domeflats, domefile, [1225,50,1265,286], $
;     display=doplot, psfile='qaplot_domeflat_01dec22.ps'
;
;   skyflatfile = 'n3.skyflat.fits'
;   skyflats = file_search('n3006[2-8].fits')
;   arm_flatcombine, skyflats, skyflatfile, [1225,50,1265,286], $
;     display=doplot, psfile='qaplot_skyflat_01dec22.ps'
;
; Dec 23 - no sky flats, see below
;
;   biasfile = 'n4.zero.fits'
;   biases = file_search(['n400[0-2]?.fits','n40030.fits'])
;   arm_zerocombine, biases, biasfile, display=doplot, psfile='qaplot_bias_01dec23.ps'
;
;   domefile = 'n4.domeflat.fits'
;   domeflats = file_search(['n4003[1-9].fits','n400[4-5]?.fits','n40060.fits'])
;   arm_flatcombine, domeflats, domefile, [1225,50,1265,286], $
;     display=doplot, psfile='qaplot_domeflat_01dec23.ps'
;   
; Dec 24 - no sky flats, see below
;
;   biasfile = 'n5.zero.fits'
;   biases = file_search(['n500[0-2]?.fits','n50030.fits'])
;   arm_zerocombine, biases, biasfile, display=doplot, psfile='qaplot_bias_01dec24.ps'
;
;   domefile = 'n5.domeflat.fits'
;   domeflats = file_search(['n5003[1-9].fits','n500[4-5]?.fits','n50060.fits'])
;   arm_flatcombine, domeflats, domefile, [1225,50,1265,286], $
;     display=doplot, psfile='qaplot_domeflat_01dec24.ps'
;
; no sky flats observed on Dec 23 or Dec 24, so make a master sky flat
; from the preceeding three nights
;
;   skyflatfile = 'skyflat_01dec.fits'
;   skyflats = file_search(['n1006[3-8].fits','n2006[6-9].fits',$
;     'n2007[0-1].fits','n3006[2-8].fits'])
;   arm_flatcombine, skyflats, skyflatfile, [1225,50,1265,286], $
;     display=doplot, psfile='qaplot_skyflat_01dec.ps'

; correct some headers

;   h = headfits('n10074.fits')
;   sxaddpar, h, 'OBJECT', 'ngc1097-20arcsec'
;   modfits, 'n10074.fits', 0, h
;
;   h = headfits('n10075.fits')
;   sxaddpar, h, 'OBJECT', 'ngc1097-20arcsec'
;   modfits, 'n10075.fits', 0, h
    
;   h = headfits('n30074.fits')
;   sxaddpar, h, 'OBJECT', 'ngc1097-nuc'
;   modfits, 'n30074.fits', 0, h

;   h = headfits('n30075.fits')
;   sxaddpar, h, 'OBJECT', 'ngc1097-nuc'
;   modfits, 'n30075.fits', 0, h
    
;   h = headfits('n50075.fits')
;   sxaddpar, h, 'OBJECT', 'ngc3351-nuc'
;   modfits, 'n50075.fits', 0, h
    
; update the headers
    
;   headfile = 'header_keywords_01dec.dat'
;   iheader_keywords, headfile, /update, silent=silent

; ---------------------------------------------------------------------------
; carry out all the initial reductions; no arc lamps were observed on
; Dec 23, so use the lamps from Dec 22
; ---------------------------------------------------------------------------
    
    iall, paramfile, procfile=procfile, tracefile=tracefile, doplot=doplot, rednight=[0,1,2,4]
    iall, paramfile, procfile=procfile, tracefile=tracefile, doplot=doplot, rednight=3, /noarcfit

; ---------------------------------------------------------------------------    
; combine cr-splits and cosmic-ray reject; see README on n50066-7
; ---------------------------------------------------------------------------    

;   ibatch, paramfile[0], crsplitfile=crsplitfile[0], /find_pixshift, /debug, /crsplits

    iall, paramfile, crsplitfile=crsplitfile, crfile=crfile, /find_pixshift

; n20074-5/NGC1097-56arcsec
    icrcombine, 'rn2'+['0074','0075']+'.fits', /find_pixshift, $
      user_pixshift=[0.0,3.0], /wfits;, /debug
    ibatch, paramfile[1], crlist='rn20074_2.fits', /crclean

; clean up residual cosmic rays

    ibatch, paramfile[0], crlist='rn10088_2.fits', /crclean, objlim=0.5, iaxis=1
    ibatch, paramfile[1], crlist='rn20097.fits', /crclean, objlim=0.5, iaxis=1
    ibatch, paramfile[2], crlist='rn30076_2.fits', /crclean, objlim=0.5, iaxis=0
    ibatch, paramfile[2], crlist='rn30078_2.fits', /crclean, objlim=1.5, iaxis=1
    ibatch, paramfile[4], crlist='rn50067.fits', /crclean, objlim=1.0, iaxis=1

; repair LA_COSMIC pixels

    repair_la_cosmic_crpix, 'crn30076_2.fits', 'repair.n30076_2.dat', /repair

; ---------------------------------------------------------------------------    
; select sky apertures and sky subtract
; ---------------------------------------------------------------------------    

;   iskyselect, skyfile=skyfile, skyapfile=skyapfile, norder_sky=1, /skyremember

    iall, paramfile, skyapfile=skyapfile

; ---------------------------------------------------------------------------    
; make the sensitivity curves
; ---------------------------------------------------------------------------    

; generate the average sensitivity curve since no standards were
; observed on 23 or 24 December

    sensname = 'sens_01dec.fits' & senstitle='2001 December'
    readcol, stdallfile, stdlist, format='A', /silent, comment='#'
    ibatch, paramfile[0], stdlist=stdlist, sensname=sensname, /makesens, zptshift=zptshift, $
      grey=2, sensinfo=sensinfo, senstitle=senstitle, slit_width=3.0, doplot=doplot, $
      wmapname='wmaplist_01dec.txt', aperture=20.0

;   info = isensfunc(i1dnames(stdlist,aperture=20.0),grey=2,sensname=sensname,$
;     senstitle=senstitle,doplot=doplot,/write,slit_width=3.0,zptshift=zptshift)

; do not generate a telluric spectrum because standards were not taken
; every night     
    
; ---------------------------------------------------------------------------    
; flux-calibrate
; ---------------------------------------------------------------------------    
    
    sensname = 'sens_01dec.fits'
    iall, paramfile, fluxcalibfile=calibfile, calibfile=calibfile, sensname=sensname

; ---------------------------------------------------------------------------    
; miscellaneous tasks
; ---------------------------------------------------------------------------    

; ---------------------------------------------------------------------------    
; generate the web page
; ---------------------------------------------------------------------------    

    if keyword_set(web) then begin
       readcol, 'weblist_01dec.txt', weblist, format='A', comment='#', /silent
       ispec_webpage, '01dec', weblist=weblist, html_path=atlas_path(/dataweb), $
         html_only=html_only, /cleanpng
    endif

; ---------------------------------------------------------------------------    
; old code
; ---------------------------------------------------------------------------    

; fix individual frames for cosmic rays

;   ibatch, paramfile[0], crlist='rn10077.fits', /crclean, objlim=1.5 ; N1097
;   ibatch, paramfile[0], crlist='rn10088_2.fits', /crclean, objlim=0.5, sigclip=3.8 ; N2915
;   ibatch, paramfile[0], crlist='rn10094.fits', /crclean, objlim=1.0 ; N3621
;   ibatch, paramfile[1], crlist='rn20074_2.fits', /crclean, objlim=1.6 ; N1097
;   ibatch, paramfile[1], crlist='rn20081_2.fits', /crclean, objlim=1.2 ; N1316
;   ibatch, paramfile[1], crlist='rn20083_2.fits', /crclean, objlim=1.2 ; N1566
;   ibatch, paramfile[1], crlist='rn20093_2.fits', /crclean, objlim=1.2, iaxis=1 ; N3621
;   ibatch, paramfile[1], crlist='rn20097.fits', /crclean, objlim=0.9, sigfrac=0.3, sigclip=3.5, iaxis=1 ; N3521
;
;   ibatch, paramfile[2], crlist='rn30076_2.fits', /crclean, objlim=0.8, sigclip=4.0 ; N7793
;   ibatch, paramfile[2], crlist='rn30078_2.fits', /crclean, objlim=1.0, iaxis=1 ; N7793
;   ibatch, paramfile[2], crlist='rn30084_2.fits', /crclean, objlim=1.5 ; N1404
;   ibatch, paramfile[2], crlist='rn30087_2.fits', /crclean, objlim=1.5 ; N1404
;   ibatch, paramfile[2], crlist='rn30095_2.fits', /crclean, objlim=1.2, iaxis=1 ; N1705
;   ibatch, paramfile[2], crlist='rn30097_2.fits', /crclean, objlim=0.8, sigclip=4.0 ; N3521
;   
;   rcube = rd2dspec('rn30084_2.fits') ; i do not believe that it's a cosmic ray jm02dec3uofa
;   ccube = rd2dspec('crn30084_2.fits')
;   ccube.image[717:724,59:67] = rcube.image[717:724,59:67]
;   ccube.mask[717:724,59:67] = fix(0)
;   wrt2dspec, ccube.fname, ccube.image, ccube.sigmamap, ccube.mask, *ccube.header
;   icleanup, rcube & icleanup, ccube
;
;   ibatch, paramfile[3], crlist='rn50065.fits', /crclean, objlim=0.8, sigclip=4.0 ; N1291
;   ibatch, paramfile[3], crlist='rn50066.fits', /crclean, objlim=1.2, sigclip=4.0 ; N1482
;   ibatch, paramfile[3], crlist='rn50067.fits', /crclean, objlim=1.0, sigclip=4.0, iaxis=1 ; N1482
;   ibatch, paramfile[3], crlist='rn50073_2.fits', /crclean, objlim=0.8, sigclip=3.8 ; N3351
;   ibatch, paramfile[3], crlist='rn50076.fits', /crclean, objlim=1.0 ; N3521

return
end    
