pro apr05_script, doplot=doplot, preproc=preproc, zptshift=zptshift, web=web
; jm06jan05uofa
    
    dateroot = '05apr'
    dates = dateroot+['04','05','06','07','08']

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

    if keyword_set(preproc) then begin
    
; verify header fidelity; 

       iheader_check, root='a.'

; clean up some headers

       flist = 'a.'+[$
         '0408',$
         '0409',$
         '0410',$
         '0411',$
         '0417',$
         '0418',$
         '0419',$
         '0420',$
         '0422',$
         '0423',$
         '0424',$
         '0425',$
         '0722',$
         '0723',$
         '0724']+'.fits'

       object = [$
         'NGC2976 E drift56',$
         'NGC2976 E drift56',$
         'NGC2976 W drift56',$
         'NGC2976 W drift56',$
         'NGC3034 S drift56',$
         'NGC3034 S drift56',$
         'NGC3034 N drift56',$
         'NGC3034 N drift56',$
         'NGC4321 SE drift56',$
         'NGC4321 SE drift56',$
         'NGC4321 NW drift56',$
         'NGC4321 NW drift56',$
         'NGC6503 drift',$
         'NGC6503 drift',$
         'NGC6503 drift']
       niceprint, flist, object
    
       for iobj = 0L, n_elements(flist)-1L do begin
          h = headfits(flist[iobj])
          sxaddpar, h, 'OBJECT', object[iobj]
          djs_modfits, flist[iobj], 0, h
       endfor
       
; average the bias frames and dome flats from each night

       root = 'a.'+['04','05','06','07','08']
       dateroot = '05apr'+['04','05','06','07','08']
       biasname = root+'01.fits'
       domename = root+'02.fits'
       biasqaplot = 'qaplot_bias_'+dateroot+'.ps'
       domeqaplot = 'qaplot_domeflat_'+dateroot+'.ps'

       overscan = [1207,0,1218,153]
       
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
          
       endfor

; combine sky flats from all the nights

       skyflats = file_search('a.0[4,5,8]0[3-5].fits')
       arm_flatcombine, skyflats, 'skyflat_05apr.fits', overscan, nhigh=1, $
         nlow=1, psfile='qaplot_skyflat_05apr.ps'

; generate the parameter file lists
       
;      ispec_makelists, 'a.'+['04','05','06','07','08'], '05apr'+['04','05','06','07','08'], $
;        /batchfile, /gzip, lampname='HeAr', extfile='kpnoextinct.dat', gain='2.1', $
;        rdnoise='5.5', trim='5 1199 0 130', overscan='1207 1218 0 130', pscale='1.66666', $
;        minwave_guess='3625', minwave_out='3625', dwave='2.75', badpixfile='badpix.dat', $
;        tracename='', biasfile='a.'+['04','05','06','07','08']+'01.fits', $
;        domefile='a.'+['04','05','06','07','08']+'02.fits', skyfile='skyflat_05apr.fits', $
;        sensname='sens_4.5_'+'05apr'+['04','05','06','07','08']+'.fits'

; update the headers with relevant keywords
       
       headfile = 'header_keywords_05apr.dat'

;      imake_keywordfile, ['SCANLEN','PA','APERTURE'], comments=['scan length [arcsec]',$
;        'slit position angle [degrees]','slit aperture [arcsec]'], values=[0.0,90.0,2.5], $
;        root='a.0???', headfile=headfile
       
       iheader_keywords, headfile, update=1, silent=silent

    endif

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

    repair_la_cosmic_crpix, 'cra.0707_2.fits', 'repair.0707_2.dat', /repair

; NGC2976 W drift56
    icrcombine, 'ra.'+['0410','0411']+'.fits', /find_pixshift, $
      user_pixshift=[0.0,-1.0], /wfits;, /debug
    ibatch, paramfile[0], crlist='ra.0410_2.fits', /crclean

; NGC2403 SE3 drift56
    icrcombine, 'ra.'+['0513','0514']+'.fits', /find_pixshift, $
      user_pixshift=[0.0,-1.0], /wfits;, /debug
    ibatch, paramfile[1], crlist='ra.0513_2.fits', /crclean

; NGC2403 NW2 drift56
    icrcombine, 'ra.'+['0517','0518']+'.fits', /find_pixshift, $
      user_pixshift=[0.0,0.0], /wfits;, /debug
    ibatch, paramfile[1], crlist='ra.0517_2.fits', /crclean
    
; NGC3938 NW drift56
    icrcombine, 'ra.'+['0616','0617']+'.fits', /find_pixshift, $
      user_pixshift=[0.0,0.0], /wfits;, /debug
    ibatch, paramfile[2], crlist='ra.0616_2.fits', /crclean
    
; NGC4242 - 11HUGS
    icrcombine, 'ra.'+['0820','0821']+'.fits', /find_pixshift, $
      user_pixshift=[0.0,-3.0], /wfits;, /debug
    ibatch, paramfile[4], crlist='ra.0820_2.fits', /crclean

; ---------------------------------------------------------------------------    
; SINGS: stitch; select sky apertures and sky-subtract
; ---------------------------------------------------------------------------    

    sings2d_stitch, 'cra.'+['0408_2','0410_2']+'.fits', outname='wcra.0408_2_stitch.fits', $
      object='NGC2976 drift56', user_pixshift=[109.0], wmapname='wmap_ra.0414.idlsave', $
      refindx=0L, minwave=3612.0, dwave=2.75, scalefactor=[1.0,1.025], debug=1, wfits=1

    sings2d_stitch, 'cra.'+['0417_2','0419_2']+'.fits', outname='wcra.0417_2_stitch.fits', $
      object='NGC3034 drift56', user_pixshift=[102.0], wmapname='wmap_ra.0421.idlsave', $
      refindx=0L, minwave=3612.0, dwave=2.75, scalefactor=[1.0,1.07], debug=1, wfits=1

    sings2d_stitch, 'cra.'+['0422_2','0424_2']+'.fits', outname='wcra.0422_2_stitch.fits', $
      object='NGC4321 drift56', user_pixshift=[108.0], wmapname='wmap_ra.0426.idlsave', $
      refindx=0L, minwave=3612.0, dwave=2.75, scalefactor=[1.0,1.01], debug=1, wfits=1

    sings2d_stitch, 'cra.'+['0432_2','0434_2','0436_2','0438_2']+'.fits', $
      outname='wcra.0432_2_stitch.fits', object='NGC5194 drift56', scalefactor=[1.0,1.0,1.025,1.0], $
      user_pixshift=[110,108.0,114.0], wmapname='wmap_ra.0440.idlsave', $
      refindx=1L, minwave=3612.0, dwave=2.75, debug=1, wfits=1

    sings2d_stitch, 'cra.'+['0513_2','0511_2','0509_2','0515_2','0517_2','0519_2']+'.fits', $
      outname='wcra.0513_2_stitch.fits', object='NGC2403 drift56', scalefactor=[1.0,1.01,1.0,1.06,0.98,0.96], $
      user_pixshift=[85.0,92.0,85.0,90.0,85.0], wmapname='wmap_ra.0508.idlsave', $
      refindx=2L, minwave=3612.0, dwave=2.75, debug=1, wfits=1
    
    sings2d_stitch, 'cra.'+['0529_2','0527_2','0531_2','0533_2']+'.fits', $
      outname='wcra.0527_2_stitch.fits', object='NGC4236 drift56', scalefactor=[1.0,0.97,1.03,1.0], $
      user_pixshift=[66.0,109.0,114.0], wmapname='wmap_ra.0537.idlsave', $
      refindx=1L, minwave=3612.0, dwave=2.75, debug=1, wfits=1

    sings2d_stitch, 'cra.'+['0606_3','0609_4']+'.fits', outname='wcra.0603_3_stitch.fits', $
      object='NGC2841 drift56', user_pixshift=[108.0], wmapname='wmap_ra.0613.idlsave', $
      refindx=0L, minwave=3612.0, dwave=2.75, scalefactor=[1.45,1.0], debug=1, wfits=1

    sings2d_stitch, 'cra.'+['0614_2','0616_2']+'.fits', outname='wcra.0614_2_stitch.fits', $
      object='NGC3938 drift56', user_pixshift=[109.0], wmapname='wmap_ra.0618.idlsave', $
      refindx=0L, minwave=3612.0, dwave=2.75, scalefactor=[1.06,1.0], debug=1, wfits=1

    sings2d_stitch, 'cra.'+['0623_2','0621_2','0625_2']+'.fits', outname='wcra.0621_2_stitch.fits', $
      object='NGC4736 drift56', user_pixshift=[110.0,110.0], wmapname='wmap_ra.0627.idlsave', $
      refindx=1L, minwave=3612.0, dwave=2.75, scalefactor=[1.0,1.0,1.035], debug=1, wfits=1

; note the use of a different wavelength map; using 0638 changes the
; dimensions of the calibrated images in SINGS2D_STITCH    
    sings2d_stitch, 'cra.'+['0634_2','0630_3','0636_2']+'.fits', outname='wcra.0634_2_stitch.fits', $
      object='NGC5055 drift56', user_pixshift=[115.0,115.0], wmapname='wmap_ra.0627.idlsave', $
      refindx=1L, minwave=3612.0, dwave=2.75, scalefactor=[1.0,1.0,0.7], debug=1, wfits=1
    
    sings2d_stitch, 'cra.'+['0710_2','0712_2']+'.fits', outname='wcra.0710_2_stitch.fits', $
      object='NGC4125 drift56', user_pixshift=[108.0], wmapname='wmap_ra.0714.idlsave', $
      refindx=0L, minwave=3612.0, dwave=2.75, scalefactor=[1.0,1.09], debug=1, wfits=1

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
    senstitle = '2005 April '+['4','5','6','7','8']+' (2.5" Slit)'
    iall, paramfile, stdfile=stdfile2_5, senstitle=senstitle, sensinfo=sensinfo, $
      slit_width=2.5, grey=2, doplot=doplot, sensname=senslist

;   readcol, stdfile2_5[i], stdlist, format='A', /silent, comment='#'
;   info = isensfunc(i1dnames(stdlist,aperture=21.83),grey=2,sensname=senslist[i],$
;     senstitle=senstitle[i],doplot=doplot,/write,slit_width=4.5,zptshift=0.0)

; generate sensitivity curves for each night (4.5")

    senslist = 'sens_4.5_'+dates+'.fits'
    senstitle = '2005 April '+['4','5','6','7','8']+' (4.5" Slit)'
    iall, paramfile, stdfile=stdfile4_5, senstitle=senstitle, sensinfo=sensinfo, $
      slit_width=4.5, grey=2, doplot=doplot

;   readcol, stdfile4_5[i], stdlist, format='A', /silent, comment='#'
;   info = isensfunc(i1dnames(stdlist,aperture=21.83),grey=2,sensname=senslist[i],$
;     senstitle=senstitle[i],doplot=doplot,/write,slit_width=4.5,zptshift=0.0)

; make the combined sensitivity functions
    
    sensname = 'sens_2.5_'+dateroot+'.fits' & senstitle='2005 April (2.5" Slit)'
    readcol, stdallfile2_5, stdlist, format='A', /silent, comment='#'
    ibatch, paramfile[0], stdlist=stdlist, sensname=sensname, /makesens, zptshift=0.0, $
      grey=2, sensinfo=sensinfo, senstitle=senstitle, slit_width=2.5, doplot=doplot, $
      wmapname='wmaplist_'+dateroot+'.txt'
    
;   info = isensfunc(i1dnames(stdlist,aperture=21.83),grey=2,sensname=sensname,$
;     senstitle=senstitle,doplot=doplot,/write,slit_width=2.5,zptshift=0.0)

    sensname = 'sens_4.5_'+dateroot+'.fits' & senstitle='2005 April (4.5" Slit)'
    readcol, stdallfile4_5, stdlist, format='A', /silent, comment='#'
    ibatch, paramfile[0], stdlist=stdlist, sensname=sensname, /makesens, zptshift=zptshift, $
      grey=2, sensinfo=sensinfo, senstitle=senstitle, slit_width=4.5, doplot=doplot, $
      wmapname='wmaplist_'+dateroot+'.txt'

;   info = isensfunc('w'+i1dnames(stdlist,aperture=21.5),grey=2,sensname=sensname,$
;     senstitle=senstitle,doplot=1,/write,slit_width=4.5,zptshift=zptshift)

; ---------------------------------------------------------------------------    
; generate the telluric spectra
; ---------------------------------------------------------------------------    

    for i = 0L, n_elements(paramfile)-1L do begin

       readcol, tellfile[i], stdlist, format='A', /silent, comment='#'
       case i of
          1L: aperture = 21.50
          4L: aperture = 22.33
          else: aperture = 21.83
       endcase

       hotlist = i1dnames(stdlist,aperture=aperture)
       iconstruct_telluric, hotlist, tellmethod=2L, contmethod=1L, method=3L, $
         npoly=2, psname='qaplot_'+repstr(tellfits[i],'.fits','.ps'), tellfits=tellfits[i], $
         doplot=doplot, debug=0, /write

    endfor
    
; ---------------------------------------------------------------------------    
; flux-calibrate
; ---------------------------------------------------------------------------    

; all the observations from 05apr06 were stitched, so skip that night    
    
    sensname = 'sens_4.5_05apr.fits'
    iall, paramfile, fluxcalibfile=calibfile, calibfile=calibfile, $
      sensname=sensname, rednight=[0,1,3,4]; , tellfits=tellfits

; the stitched spectra have been wavelength calibrated already, so
; just flux calibrate    
    
    readcol, 'skylist_stitch.txt', caliblist, format='A', comment='#', /silent
    ibatch, paramfile[0], caliblist='s'+caliblist, $ ; tellfits=tellfits[0], $
      sensname=sensname, /fluxcalibrate

; ---------------------------------------------------------------------------    
; miscellaneous tasks
; ---------------------------------------------------------------------------    

; compare the sensitivity curves

    label = [['April 4','Clear, Windy'],['April 5','Clear'],['April 6','Partly Cloudy'],$
      ['April 7','Partly Cloudy, Windy'],['April 8','Clear'],['2.5" Slit',''] ]
    senslist = ['sens_4.5_'+dates+'.fits','sens_2.5_'+dateroot+'.fits']
    meansens = 'sens_4.5_'+dateroot+'.fits'
    isenscompare, senslist, meansens=meansens, label=label, title='2005 April', $
      psname='qaplot_sens_compare_'+dateroot+'.ps', /postscript

    loglist = ['qalog_sens_4.5_'+dates+'.log','qalog_sens_'+['2.5','4.5']+'_'+dateroot+'.log']
    label = [['April 4','Clear, Windy'],['April 5','Clear'],['April 6','Partly Cloudy'],$
      ['April 7','Partly Cloudy, Windy'],['April 8','Clear'],['2.5" Slit',''],['4.5" Slit',''] ]
    iqaphoto, loglist, label=label, title='2005 April', psname='qaplot_photo_'+dateroot+'.ps', /postscript
    
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
