pro may06_script, doplot=doplot, preproc=preproc, zptshift=zptshift, web=web
; jm06jul07uofa

    dateroot = '06may'
    dates = dateroot+['21','22']

    paramfile = 'ibatch_'+dates+'.txt'
    procfile = 'objlist_'+dates+'.txt'
    skyfile = 'skylist_'+dates+'.txt'
    skyapfile = 'skyaplist_'+dates+'.txt'
    crsplitfile = 'crsplits_'+dates+'.txt'
    crfile = 'crlist_'+dates+'.txt'
    tracefile = 'tracelist_'+dates+'.txt'
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

    if keyword_set(preproc) then begin
    
; verify header fidelity; 

       iheader_check, root='a.'

; average the bias frames, dome flats, and sky flats from each night 

       fileroot = ['11','12']
       root = 'a.'+fileroot
       biasname = 'bias_'+dates+'.fits'
       domename = 'domeflat_'+dates+'.fits'
       skyflatname = 'skyflat_'+dates+'.fits'
       biasqaplot = 'qaplot_bias_'+dates+'.ps'
       domeqaplot = 'qaplot_domeflat_'+dates+'.ps'
       skyflatqaplot = 'qaplot_skyflat_'+dates+'.ps'

; combine sky flats from all the nights

       overscan = [1207,4,1218,141]

       for iroot = 0L, n_elements(root)-1L do begin

          biases = file_search(root[iroot]+'90[0-2]?.fits')
          arm_zerocombine, biases, biasname[iroot], nhigh=1, nlow=1, $
            psfile=biasqaplot[iroot]

          domeflats = file_search(root[iroot]+'90[5-6]?.fits')
          arm_flatcombine, domeflats, domename[iroot], overscan, nhigh=1, nlow=1, $
            psfile=domeqaplot[iroot]

; estimate the gain and the read noise of each night

          flatindx = fix(randomu(seed,1L)*n_elements(domeflats)-1L)
          flatfile1 = domeflats[flatindx] & flatfile2 = domeflats[flatindx+1L]
          biasindx = fix(randomu(seed,1L)*n_elements(biases)-1L)
          biasfile1 = biases[biasindx] & biasfile2 = biases[biasindx+1L]
          
;         gain_rdnoise, flatfile1[0], flatfile2[0], biasfile1[0], biasfile2[0], $
;           xskip=xskip, yave=yave, gain=gain, rdnoise=rdnoise
;         splog, 'Night '+string(iroot+1L,format='(I0)')+': Gain = '+$
;           strtrim(string(gain,format='(F12.1)'),2)+$
;           ', readnoise = '+strtrim(string(rdnoise,format='(F12.1)'),2)

; master sky flat          

          if (iroot eq 0L) then $
            skyflats = file_search(root[iroot]+'0[1-9].fits') else $
            skyflats = file_search(root[iroot]+'0[1-7].fits')
            
          arm_flatcombine, skyflats, skyflatname[iroot], overscan, nhigh=1, $
            nlow=1, psfile=skyflatqaplot[iroot]

       endfor

; generate the parameter file lists

;      ispec_makelists, 'a.'+fileroot, dates, /batchfile, gzip=0, lampname='HeAr', $
;        extfile='kpnoextinct.dat', gain='2.1', rdnoise='5.5', trim='5 1197 4 141', $
;        overscan='1207 1218 4 141', pscale='1.66666', minwave_guess='3585', $
;        minwave_out='3590', dwave='2.75', badpixfile='badpix.dat', biasfile=biasname, $
;        domefile=domename, skyfile=skyflatname, sensname='sens_4.5_'+dates+'.fits', /all
 
; update the headers with relevant keywords
       
       headfile = 'header_keywords_'+dateroot+'.dat'
 
;      imake_keywordfile, ['SCANLEN','POSANGLE','APERTURE'], comments=['scan length [arcsec]',$
;        'slit position angle [degrees]','slit aperture [arcsec]'], values=[0.0,90.0,2.5], $
;        root='a.1???', headfile=headfile

       iheader_keywords, headfile, update=1, silent=silent

    endif

; ---------------------------------------------------------------------------
; carry out all the initial reductions
; ---------------------------------------------------------------------------
    
    iall, paramfile, procfile=procfile, tracefile=tracefile, doplot=doplot

; ---------------------------------------------------------------------------    
; combine cr-splits and cosmic-ray reject
; ---------------------------------------------------------------------------    

;   ibatch, paramfile[0], crsplitfile=crsplitfile[0], /find_pixshift, /debug, /crsplits
    iall, paramfile, crsplitfile=crsplitfile, crfile=crfile, /find_pixshift

; ---------------------------------------------------------------------------    
; SINGS: stitch; select sky apertures and sky-subtract
; ---------------------------------------------------------------------------    

    sings2d_stitch, 'cra.11'+['41','39','43','45']+'_2.fits', outname='wcra.1139_2_stitch.fits', $
      object='NGC6946 drift56', scalefactor=[1.15,1.15,1.0,1.0], $
      user_pixshift=[109.0,108.0,109.0], wmapname='wmap_ra.1148.idlsave', $
      refindx=1L, minwave=3620.0, dwave=2.75, debug=1, wfits=1

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

    senslist = 'sens_2.5_'+dates+'.fits'
    senstitle = '2006 May '+['21','22']+' (2.5" Slit)'
    iall, paramfile, stdfile=stdfile2_5, senstitle=senstitle, sensinfo=sensinfo, $
      slit_width=2.5, grey=2, doplot=doplot, sensname=senslist

;   readcol, stdfile2_5[i], stdlist, format='A', /silent, comment='#'
;   info = isensfunc(i1dnames(stdlist,aperture=23),grey=2,sensname=senslist[i],$
;     senstitle=senstitle[i],doplot=1,/write,slit_width=4.5,zptshift=0.0)

; generate sensitivity curves for each night (4.5")

    senslist = 'sens_4.5_'+dates+'.fits'
    senstitle = '2006 May '+['21','22']+' (4.5" Slit)'
    iall, paramfile, stdfile=stdfile4_5, senstitle=senstitle, sensinfo=sensinfo, $
      slit_width=4.5, grey=2, doplot=doplot

;   readcol, stdfile4_5[i], stdlist, format='A', /silent, comment='#'
;   info = isensfunc(i1dnames(stdlist,aperture=23),grey=2,sensname=senslist[i],$
;     senstitle=senstitle[i],doplot=1,/write,slit_width=4.5,zptshift=0.0)

; make the combined sensitivity functions
    
    sensname = 'sens_2.5_'+dateroot+'.fits' & senstitle='2006 May  (2.5" Slit)'
    readcol, stdallfile2_5, stdlist, format='A', /silent, comment='#'
    ibatch, paramfile[0], stdlist=stdlist, sensname=sensname, /makesens, zptshift=0.0, $
      grey=2, sensinfo=sensinfo, senstitle=senstitle, slit_width=2.5, doplot=doplot, $
      wmapname='wmaplist_'+dateroot+'.txt'
    
;   info = isensfunc(i1dnames(stdlist,aperture=23),grey=2,sensname=sensname,$
;     senstitle=senstitle,doplot=1,/write,slit_width=2.5,zptshift=0.0)

    sensname = 'sens_4.5_'+dateroot+'.fits' & senstitle='2006 May (4.5" Slit)'
    readcol, stdallfile4_5, stdlist, format='A', /silent, comment='#'
    ibatch, paramfile[0], stdlist=stdlist, sensname=sensname, /makesens, zptshift=zptshift, $
      grey=2, sensinfo=sensinfo, senstitle=senstitle, slit_width=4.5, doplot=doplot, $
      wmapname='wmaplist_'+dateroot+'.txt'

;   info = isensfunc(i1dnames(stdlist,aperture=23),grey=2,sensname=sensname,$
;     senstitle=senstitle,doplot=1,/write,slit_width=4.5,zptshift=zptshift)

; ---------------------------------------------------------------------------    
; generate the telluric spectra
; ---------------------------------------------------------------------------    

    for i = 0L, n_elements(paramfile)-1L do begin
       readcol, tellfile[i], stdlist, format='A', /silent, comment='#'
       hotlist = i1dnames(stdlist,aperture=23.0)
       iconstruct_telluric, hotlist, tellmethod=2L, contmethod=1L, method=3L, $
         npoly=2, psname='qaplot_'+repstr(tellfits[i],'.fits','.ps'), tellfits=tellfits[i], $
         doplot=doplot, debug=0, /write
    endfor
    
; ---------------------------------------------------------------------------    
; flux-calibrate and rectify the HII-region spectra (which have not
; been sky-subtracted) 
; ---------------------------------------------------------------------------    

    sensname = 'sens_4.5_'+dateroot+'.fits'
    iall, paramfile, calibfile=calibfile, sensname=sensname;, tellfits=tellfits

; the stitched spectrum has been wavelength calibrated already, so
; just flux calibrate
    
    sensname = 'sens_4.5_'+dateroot+'.fits'
    readcol, 'skylist_stitch.txt', caliblist, format='A', comment='#', /silent
    ibatch, paramfile[0], caliblist='s'+caliblist, $ ; tellfits=tellfits[0], $
      sensname=sensname, /fluxcalibrate

; ---------------------------------------------------------------------------    
; miscellaneous tasks
; ---------------------------------------------------------------------------    

; compare the sensitivity curves

    label = [['May 21','Mostly clear, windy'],['May 22','Clear!'],['2.5" Slit','']]
    senslist = ['sens_4.5_'+dates+'.fits','sens_2.5_'+dateroot+'.fits']
    meansens = 'sens_4.5_'+dateroot+'.fits'
    isenscompare, senslist, meansens=meansens, label=label, title='2006 May', $
      psname='qaplot_sens_compare_'+dateroot+'.ps', /postscript

    loglist = ['qalog_sens_4.5_'+dates+'.log','qalog_sens_'+['2.5','4.5']+'_'+dateroot+'.log']
    label = [['May 21','Mostly clear, windy'],['May 22','Clear!'],['2.5" Slit',''],['4.5" Slit','']]
    iqaphoto, loglist, label=label, title='2006 May', psname='qaplot_photo_'+dateroot+'.ps', /postscript
    
; ---------------------------------------------------------------------------    
; generate the web page
; ---------------------------------------------------------------------------    

    if keyword_set(web) then begin
       readcol, 'weblist_'+dateroot+'.txt', weblist, format='A', comment='#', /silent
       ispec_webpage, dateroot, weblist=weblist, html_path=atlas_path(/dataweb), $
         html_only=html_only, /cleanpng
    endif

return 
end   
