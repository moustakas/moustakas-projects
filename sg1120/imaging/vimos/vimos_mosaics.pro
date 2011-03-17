pro vimos_mosaics, dec03=dec03, feb06=feb06, sextractor=sextractor, scamp=scamp, $
  swarp=swarp, test=test, jpeg=jpeg, use_both_epochs=use_both_epochs, band=band
; jm07jun - generate mosaics from the VIMOS 2003 and 2006 data; by
;           setting USE_BOTH_EPOCHS=1 (recommended) you can run SCAMP and
;           SWARP on both epoch simultaneously 
; jm08aug - major update; astrometry for both epochs is now good
;           to 0.04" (rms, internal) and 0.3" (rms, reference); XML
;           pages are also now built and the configuration structures
;           are saved

; see README    
    
    if (not keyword_set(use_both_epochs)) then $
      if ((not keyword_set(dec03)) and (not keyword_set(feb06))) or $
      (keyword_set(dec03) and keyword_set(feb06)) then begin
       splog, 'Must specify either DEC03, FEB06, or USE_BOTH_EPOCHS keyword!'
       return
    endif
    
    if (keyword_set(use_both_epochs) and keyword_set(dec03)) or $
      (keyword_set(use_both_epochs) and keyword_set(feb06)) then begin
       splog, 'You cannot simultaneously use USE_BOTH_EPOCHS and DEC03 or FEB06 keywords!'
       return
    endif
    
    sexpath = sg1120_path(/sex)
    dec03_datapath = vimos_path(dec03=1,feb06=0)+'sg1120/'
    feb06_datapath = vimos_path(dec03=0,feb06=1)+'sg1120/'
    datapath = vimos_path(dec03=dec03,feb06=feb06)+'sg1120/'
    use_both_epochs_datapath = vimos_path(/mosaics)+'data/'
    mosaicpath = vimos_path(/mosaics)

    suffix = ''
    if keyword_set(dec03) then suffix = '_2003'
    if keyword_set(feb06) then suffix = '_2006'

    if (not keyword_set(band)) then band = '[B,V,R]'
    bandpass = strsplit(repstr(repstr(band,'[',''),']',''),',',/extract)
    nband = n_elements(bandpass)
    gain = string([1.73,1.86,1.95,1.80],format='(F4.2)') ; [e/ADU]

;   mosaicpath = vimos_path(/mosaics)+'test/'
;   datapath = mosaicpath

    stiffconfig = sexpath+'default.stiff'
    scampconfig = sexpath+'default.scamp'
    swarpconfig = sexpath+'default.swarp'
    sexconfig = sexpath+'default.sex'
    sexparam = sexpath+'sg1120.sex.param'
    sexconv = sexpath+'default.conv'
    sexnnw = sexpath+'default.nnw'

    http = 'http://sdss.physics.nyu.edu/ioannis/research/sg1120/sex/'
    xslsex = http+'sex.xsl'
    xslscamp = http+'scamp.xsl'
    xslswarp = http+'swarp.xsl'

; ---------------------------------------------------------------------------    
; generate SE catalogs
    if keyword_set(sextractor) then begin
       for ib = 0, nband-1 do begin
          if keyword_set(use_both_epochs) then begin
             splog, 'SEXTRACTOR keyword incompatible '+$
               'with USE_BOTH_EPOCHS keyword'
             return
          endif

          imagelist = file_search(datapath+'ra.sg1120*_'+$
            bandpass[ib]+'.fits',count=nimage)
          weightlist = repstr(imagelist,'.fits','.weight.fits')
          rmslist = repstr(imagelist,'.fits','.rms.fits')
          flaglist = repstr(imagelist,'.fits','.flag.fits')

          catlist = repstr(imagelist,'.fits','.cat')
          radecreglist = repstr(catlist,'.cat','.reg')
          seglist = repstr(imagelist,'.fits','.seg.fits')

          info = im_headerforage(imagelist,ext=1) ; grab header info
       
; initialize the SE configuration parameters
          config = init_sex_config(nimage)
          configfile = datapath+'sex.config'
       
          config.catalog_name = catlist
          config.weight_image = rmslist
          config.flag_image = flaglist
          config.parameters_name = sexparam
          config.filter_name = sexconv
          config.starnnw_name = sexnnw

          config.catalog_type = 'FITS_LDAC'
          config.detect_thresh = 1.5
          config.analysis_thresh = 1.5
          config.weight_type = 'MAP_RMS'
          config.weight_gain = 'N'
          config.interp_type = 'NONE'
          config.nthreads = 4

          config.seeing_fwhm = info.seeing
          config.mag_zeropoint = 28.0
;         config.mag_zeropoint = info.mag0 ; use the magnitude zeropoint in scamp, not here
       
          config.checkimage_type = 'NONE' ; SEGMENTATION
          config.checkimage_name = seglist

          mwrfits, config, configfile+'.fits', /create

; do it!       
          t0 = systime(1)
          im_sex, imagelist, config, silent=silent ; do not pass CONFIGFILE
          splog, 'Total time to generate SE catalogs = ', $
            (systime(1)-t0)/60.0, ' minutes'

; make a QAplot and write out region files       
          sexqaplot = mosaicpath+'qaplot_'+bandpass[ib]+'_sex'+suffix+'.ps'
          dfpsplot, sexqaplot
          for ii = 0L, nimage-1L do begin
             for jj = 0L, 3L do begin
                cat1 = mrdfits(catlist[ii],2*jj+2,/silent)
                if (jj eq 0L) then cat = cat1 else cat = [cat,cat1]
             endfor
             xr = [10,24] & yrange1 = [-0.02,1.0] & yrange2 = [0.0,4.9]
             plot, [0], [0], /nodata, position=[0.14,0.55,0.96,0.93], xsty=1, ysty=1, $
               xrange=xr, yrange=yrange1, xtitle='', xtickname=replicate(' ',10), $
               ytitle='Class Star', charsize=1.8, xthick=5.0, $
               ythick=5.0, charthick=5.0, title=file_basename(repstr(imagelist[ii],'.fits',''))
             djs_oplot, cat1.mag_auto, cat1.class_star, psym=symcat(16)
             djs_plot, [0], [0], /nodata, /noerase, position=[0.14,0.12,0.96,0.55], $
               xsty=1, ysty=1, xrange=xr, yrange=yrange2, xtitle='Instrumental Magnitude', $
               ytitle='r_{eff} (arcsec)', charsize=1.8, xthick=5.0, $
               ythick=5.0, charthick=5.0
             djs_oplot, cat1.mag_auto, cat1.flux_radius*cat1.awin_image*0.188, $
               psym=symcat(16), symsize=0.8

             splog, 'Witing DS9 region file '+file_basename(radecreglist[ii])
             write_ds9_regionfile, cat.xwin_world, cat.ywin_world, $
               filename=radecreglist[ii], color='red'
          endfor
          dfpsclose
       endfor
    endif

; ---------------------------------------------------------------------------
; SCAMP
    if keyword_set(scamp) then begin
       if keyword_set(use_both_epochs) then begin
          splog, 'Delete and remake all soft links in '+use_both_epochs_datapath+' ?'
          cc = get_kbrd(1)
          if strmatch(cc,'y',/fold) then begin
             print
             spawn, '/bin/rm -f '+use_both_epochs_datapath+'*'
; 2003
             imagelist = file_search(dec03_datapath+'ra.sg1120_*_'+band+'.fits',count=nimage)
             weightlist = repstr(imagelist,'.fits','.weight.fits')
             catlist = repstr(imagelist,'.fits','.cat')
             for ii = 0L, nimage-1L do spawn, 'ln -sfv '+imagelist[ii]+' '+$
               use_both_epochs_datapath+'vimos.2003.'+file_basename(imagelist[ii])
             for ii = 0L, nimage-1L do spawn, 'ln -sfv '+weightlist[ii]+' '+$
               use_both_epochs_datapath+'vimos.2003.'+file_basename(weightlist[ii])
             for ii = 0L, nimage-1L do spawn, 'ln -sfv '+catlist[ii]+' '+$
               use_both_epochs_datapath+'vimos.2003.'+file_basename(catlist[ii])
             print
; 2006
             imagelist = file_search(feb06_datapath+'ra.sg1120_*_'+band+'.fits',count=nimage)
             weightlist = repstr(imagelist,'.fits','.weight.fits')
             catlist = repstr(imagelist,'.fits','.cat')
             for ii = 0L, nimage-1L do spawn, 'ln -sfv '+imagelist[ii]+' '+$
               use_both_epochs_datapath+'vimos.2006.'+file_basename(imagelist[ii])
             for ii = 0L, nimage-1L do spawn, 'ln -sfv '+weightlist[ii]+' '+$
               use_both_epochs_datapath+'vimos.2006.'+file_basename(weightlist[ii])
             for ii = 0L, nimage-1L do spawn, 'ln -sfv '+catlist[ii]+' '+$
               use_both_epochs_datapath+'vimos.2006.'+file_basename(catlist[ii])
             print
          endif

          catlist = file_search(use_both_epochs_datapath+$
            'vimos.200?.ra*'+band+'.cat',count=nimage)

          if (nimage eq 0L) then begin
             splog, 'No catalogs found in '+use_both_epochs_datapath
             return
          endif

       endif else begin

          catlist = file_search(datapath+'ra.sg1120*_'+$
            band+'.cat',count=nimage)
          if (nimage eq 0L) then begin
             splog, 'No catalogs found in '+datapath
             return
          endif

       endelse

; initialize the scamp configuration parameters
       config = init_scamp_config()
       configfile = mosaicpath+'scamp'+suffix+'.config'

; this uses the USNO catalog built for me by Dustin       
       config.astref_catalog = 'FILE'
       config.astrefcat_name = sexpath+'sg1120_sdss_dr7_refcat.cat'
;      config.astrefcat_name = sexpath+'sg1120_usnob_refcat.cat'
       config.astrefcent_keys = 'RA,DEC'
       config.astreferr_keys = 'ERR_A,ERR_B'
       config.astrefmag_key = 'MAG'

; boost the weights
       config.astref_weight = '100.0'

;      config.astref_catalog = 'NOMAD-1' & config.astref_band = 'B'
;      config.astref_catalog = '2MASS'   & config.astref_band = 'Ks' 
;      config.astref_catalog = 'USNO-B1' & config.astref_band = 'Bj' 
       
       config.crossid_radius = 1.5 ; this is important

; the OBJECT names are different between 2003 & 2006, so use this
; header tag to force scamp to define a different *astrometric*
; instrument between the two epochs; use the FILTER keyword to define
; the three (BVR) photometric instruments
       
       config.astrinstru_key = 'FILTER,OBJECT'
       config.photinstru_key = 'FILTER'
       config.magzero_key = 'MAG0' ; 'MAGZERO'
       config.save_refcatalog = 'N'
       config.refout_catpath = mosaicpath
       config.mergedoutcat_type = 'NONE'

       config.checkplot_type = strjoin(['ASTR_CHI2','ASTR_INTERROR1D','ASTR_INTERROR2D',$
         'ASTR_REFERROR1D','ASTR_REFERROR2D','DISTORTION','FGROUPS','PHOT_ERROR','PHOT_ZPCORR',$
         'PHOT_ZPCORR3D'],',')
       config.checkplot_name = strjoin('vimos'+suffix+'_'+['astr_chi2','astr_interror1d',$
         'astr_interror2d','astr_referror1d','referror2d','distort','fgroups',$
         'psphot_error','phot_zpcorr','phot_zpcorr3d'],',')

       config.xml_name = mosaicpath+'vimos'+suffix+'.scamp.xml'
       config.xsl_url = xslscamp

       t0 = systime(1)
       maxiter = 2
       if keyword_set(dec03) then maxiter = 2
       if keyword_set(feb06) then maxiter = 2

; think carefully before increasing DEGREE; for example DEGREE>3 when
; making a single mosaic results in DOOM!       

       for iter = 0L, maxiter do begin 
          case iter of
             0L: begin
                config.distort_degrees = '1'
                config.mosaic_type = 'FIX_FOCALPLANE' ; 'UNCHANGED' ; 'LOOSE'
                config.pixscale_maxerr = '1.1'
                config.position_maxerr = '3.0'
                config.posangle_maxerr = '2.0'
                config.aheader_suffix = '.ahead'
             end
             1L: begin
                config.distort_degrees = '3,3'
                config.mosaic_type = 'FIX_FOCALPLANE' ; 'UNCHANGED' ; 'LOOSE'
                config.pixscale_maxerr = '1.1' ; '1.1'
                config.position_maxerr = '0.5' ; '0.5'
                config.posangle_maxerr = '2.0' ; '1.0'
                config.aheader_suffix = '.head'
             end
             else: begin
                config.distort_degrees = '5,5'
                config.mosaic_type = 'FIX_FOCALPLANE' ; 'UNCHANGED' ; 'LOOSE'
                config.pixscale_maxerr = '1.1'
                config.position_maxerr = '0.5'
                config.posangle_maxerr = '1.0'
                config.aheader_suffix = '.head'
             end
          endcase

          mwrfits, config, configfile+'.fits', /create
          im_scamp, catlist, config, configfile=configfile, silent=silent

       endfor 
       splog, 'Total time to run scamp = ', (systime(1)-t0)/60.0, ' minutes.'
    endif

; ---------------------------------------------------------------------------
; SWARP
    if keyword_set(swarp) then begin
; initialize the swarp configuration parameters
       config = init_swarp_config()
       if keyword_set(test) then begin
          config.center = '00:00:00.0,+00:00:00.0'
          config.center_type = 'ALL' 
          config.pixelscale_type = 'MANUAL' 
          config.pixel_scale = '0.205'
          config.image_size = '0'
          config.header_only = 'Y'
       endif else begin
          config.center_type = 'MANUAL' 
          config.pixelscale_type = 'MANUAL' 
          config.pixel_scale = '0.205'
          if keyword_set(dec03) then begin
             config.center = '11:20:03.0,-12:04:30.0'
             config.image_size = '5800,5100'
          endif
          if keyword_set(feb06) then begin
             config.center = '11:19:54.0,-12:02:43.0'
             config.image_size = '5800,5100'
          endif
          if keyword_set(use_both_epochs) then begin
             config.center = '11:19:58.0,-12:03:33.0'
             config.image_size = '6570,5670'
          endif
       endelse
       
       config.blank_badpixels = 'Y'
       config.interpolate = 'N'
       config.write_fileinfo = 'Y'
       config.celestial_type = 'EQUATORIAL'
       config.copy_keywords = 'OBJECT,FILTER'
       config.xml_name = mosaicpath+'vimos'+suffix+'.swarp.xml'
       config.xsl_url = xslswarp

; build the individual B-, V-, and R-band mosaics; note that using the
; RMS images as weight maps produces some funky results
       
       config.write_xml = 'Y'
       config.subtract_back = 'Y'
       config.resampling_type = 'LANCZOS3'
       config.weight_type = 'MAP_WEIGHT'

       config.combine_type = 'WEIGHTED'
       
       for ib = 0L, nband-1L do begin

          configfile = mosaicpath+'swarp_'+bandpass[ib]+suffix+'.config'

          config.imageout_name = mosaicpath+'sg1120_'+bandpass[ib]+suffix+'.fits'
          config.weightout_name = repstr(config.imageout_name,'.fits','.weight.fits')

          if keyword_set(use_both_epochs) then $
            imagelist = file_search(use_both_epochs_datapath+$
            'vimos.200?.ra.sg1120*_'+bandpass[ib]+'.fits') else $
              imagelist = file_search(datapath+'ra.sg1120*_'+bandpass[ib]+'.fits')
          
          rmslist = repstr(imagelist,'.fits','.rms.fits')
          weightlist = repstr(imagelist,'.fits','.weight.fits')

          config.weight_image = strjoin(weightlist,',') ; strjoin(rmslist,',') ; 

          mwrfits, config, configfile+'.fits', /create
          
          splog, 'Building '+strtrim(config.imageout_name,2)
          t0 = systime(1)
          im_swarp, imagelist, config, configfile=configfile, silent=silent
          splog, 'Total time to build '+config.imageout_name+' = ', $
            (systime(1)-t0)/60.0, ' minutes.'
          
          if (not keyword_set(test)) then begin
             fixme = [config.imageout_name,config.weightout_name]
             for ff = 0L, n_elements(fixme)-1L do begin
                hdr = headfits(fixme[ff])
;               sxaddpar, hdr, 'OBJECT', 'SG1120'
                sxaddpar, hdr, 'OBJECT', 'SG1120 '+bandpass[ib]+suffix
                modfits, fixme[ff], 0, hdr
             endfor
          endif

;; code to build the FLAGS mosaic; it doesn't properly use the
;; external scamp headers
;
;          config.write_xml = 'N'
;          config.subtract_back = 'N'
;          config.resampling_type = 'NEAREST'
;          config.combine_type = 'SUM'
;          config.weight_type = 'NONE'
;          config.header_suffix = '.flag.head'
;          
;          imagelist = repstr(imagelist,'.fits','.flag.fits')
;          config.imageout_name = repstr(config.imageout_name,$
;            '.fits','.flag.fits')
;          
;          splog, 'Building '+strtrim(config.imageout_name,2)
;          im_swarp, imagelist, config, silent=silent

       endfor 

; now build the chi2 image of all the filters *only* if
; USE_BOTH_EPOCHS=1 (since we are unlikely to do source detection on
; the individual 2003/2006 mosaics); could also do this separately for
; each filter

       if keyword_set(use_both_epochs) then begin
          
          config.combine_type = 'CHI2'
          config.copy_keywords = 'OBJECT'
          
          config.imageout_name = mosaicpath+'sg1120_BVR_chi2'+suffix+'.fits'
          config.weightout_name = repstr(config.imageout_name,'.fits','.weight.fits')

          if keyword_set(use_both_epochs) then $
            imagelist = file_search(use_both_epochs_datapath+$
            'vimos.200?.ra.sg1120*_'+band+'.fits') else $
              imagelist = file_search(datapath+'ra.sg1120*_'+band+'.fits')

          weightlist = repstr(imagelist,'.fits','.weight.fits')
          config.weight_image = strjoin(weightlist,',')

          configfile = mosaicpath+'swarp_BVR_chi2'+suffix+'.config'
          mwrfits, config, configfile+'.fits', /create

          splog, 'Building '+strtrim(config.imageout_name,2)
          t0 = systime(1)
          im_swarp, 'data/'+file_basename(imagelist), config, $
            configfile=configfile, silent=silent
          splog, 'Total time to build '+config.imageout_name+' = ', $
            (systime(1)-t0)/60.0, ' minutes.'

          if (not keyword_set(test)) then begin
             fixme = [config.imageout_name,config.weightout_name]
             for ff = 0L, n_elements(fixme)-1L do begin
                hdr = headfits(fixme[ff])
                sxaddpar, hdr, 'OBJECT', 'SG1120 BVR'+suffix+' chi2'
                modfits, fixme[ff], 0, hdr
             endfor
          endif

       endif
          
    endif 

; ---------------------------------------------------------------------------
; finally generate the color mosaics
; ---------------------------------------------------------------------------

    if keyword_set(jpeg) then begin
       tifffile = mosaicpath+'sg1120_BVR'+suffix+'.tiff'
       jpegfile = mosaicpath+'sg1120_BVR'+suffix+'.jpeg'
; JPEG
       hdr = headfits(mosaicpath+'sg1120_B'+suffix+'.fits')
       imsize = [sxpar(hdr,'NAXIS1'),sxpar(hdr,'NAXIS2')]
       nx = 512*(min(imsize)/512)

       splog, 'Reading '+mosaicpath+'sg1120_B'+suffix+'.fits'
       Bim = (mrdfits(mosaicpath+'sg1120_B'+suffix+'.fits',0,/silent))$
         [imsize[0]/2L-nx/2L:imsize[0]/2L+nx/2L-1L,$
         imsize[1]/2L-nx/2L:imsize[1]/2L+nx/2L-1L]
       splog, 'Reading '+mosaicpath+'sg1120_V'+suffix+'.fits'
       Vim = (mrdfits(mosaicpath+'sg1120_V'+suffix+'.fits',0,/silent))$
         [imsize[0]/2L-nx/2L:imsize[0]/2L+nx/2L-1L,$
         imsize[1]/2L-nx/2L:imsize[1]/2L+nx/2L-1L]
       splog, 'Reading '+mosaicpath+'sg1120_R'+suffix+'.fits'
       Rim = (mrdfits(mosaicpath+'sg1120_R'+suffix+'.fits',0,/silent))$
         [imsize[0]/2L-nx/2L:imsize[0]/2L+nx/2L-1L,$
         imsize[1]/2L-nx/2L:imsize[1]/2L+nx/2L-1L]

       splog, 'Writing '+jpegfile
       scales = [0.9,1.5,3.5]*1D11
;      scales = [1.0,1.8,2.5]*1D11
       nw_rgb_make, Rim, Vim, Bim, name=jpegfile, scales=scales, $
         nonlinearity=3.0, rebinfactor=2, quality=75

; TIFF
       splog, 'Writing '+tifffile
       mosaic_file = strjoin(mosaicpath+'sg1120_'+['R','V','B']+suffix+'.fits',' ')
       spawn, 'stiff -c '+stiffconfig+' '+mosaic_file+' -OUTFILE_NAME '+tifffile+$
         ' -BINNING 2 -GAMMA_FAC 0.9 -COLOUR_SAT 1.8'
    endif

; ###########################################################################
; combine the 2003+2006 mosaics into a single mosaic
; ###########################################################################

;;; ---------------------------------------------------------------------------
;;; onemosaic/sextractor
;;; ---------------------------------------------------------------------------
;;
;;    if keyword_set(onemosaic_sextractor) then begin
;;
;;       imagelist = file_search(mosaicpath+'vimos_'+band+$
;;         '_200[3,6].fits',count=nimage)
;;       weightlist = repstr(imagelist,'.fits','.weight.fits')
;;       catlist = repstr(imagelist,'.fits','.cat')
;;       radecreglist = repstr(catlist,'.cat','.reg')
;;
;;       config = init_sex_config(nimage)
;;       
;;       config.catalog_name = catlist
;;       config.weight_image = weightlist
;;       config.parameters_name = sexparam
;;       config.filter_name = sexconv
;;       config.starnnw_name = sexnnw
;;
;;       config.catalog_type = 'FITS_LDAC'
;;       config.detect_thresh = 3.0
;;       config.analysis_thresh = 3.0
;;       config.weight_type = 'MAP_WEIGHT'
;;       config.weight_gain = 'N'
;;       config.interp_type = 'NONE'
;;       config.nthreads = 4
;;       config.seeing_fwhm = 0.7
;;
;;       t0 = systime(1)
;;       im_sex, imagelist, config, silent=silent
;;       splog, 'Total time to generate SE catalogs = ', $
;;         (systime(1)-t0)/60.0, ' minutes'
;;
;;; write out some handy DS9 region files       
;;
;;       for ii = 0L, nimage-1L do begin
;;          cat = mrdfits(catlist[ii],2,/silent)
;;          splog, 'Witing DS9 region file '+file_basename(radecreglist[ii])
;;          write_ds9_regionfile, cat.xwin_world, cat.ywin_world, $
;;            filename=radecreglist[ii], color='red'
;;       endfor
;;
;;    endif
;;
;;; ---------------------------------------------------------------------------
;;; onemosaic/scamp
;;; ---------------------------------------------------------------------------
;;
;;    if keyword_set(onemosaic_scamp) then begin
;;
;;       catlist = file_search(mosaicpath+'vimos_'+band+$
;;         '_200[3,6].cat',count=ncat)
;;
;;; initialize the scamp configuration parameters
;;
;;       config = init_scamp_config()
;;       configfile = mosaicpath+'vimos'+suffix+'.scamp.config'
;;
;;       config.astref_catalog = 'FILE'
;;       config.astrefcat_name = sexpath+'sg1120_usnob_refcat.cat'
;;       config.astrefcent_keys = 'RA,DEC'
;;       config.astreferr_keys = 'ERR_A,ERR_B'
;;       config.astrefmag_key = 'MAG'
;;;      config.crossid_radius = 1.5
;;
;;       config.astrinstru_key = 'FILTER'
;;       config.photinstru_key = 'FILTER'
;;       config.save_refcatalog = 'N'
;;       config.refout_catpath = mosaicpath
;;       config.mergedoutcat_type = 'NONE'
;;
;;       config.checkplot_type = strjoin(['ASTR_CHI2','ASTR_INTERROR1D','ASTR_INTERROR2D',$
;;         'ASTR_REFERROR1D','ASTR_REFERROR2D','DISTORTION','FGROUPS','PHOT_ERROR','PHOT_ZPCORR',$
;;         'PHOT_ZPCORR3D'],',')
;;       config.checkplot_name = strjoin('vimos'+suffix+'_'+['astr_chi2','astr_interror1d',$
;;         'astr_interror2d','astr_referror1d','referror2d','distort','fgroups',$
;;         'psphot_error','phot_zpcorr','phot_zpcorr3d'],',')
;;
;;       config.xml_name = mosaicpath+'vimos'+suffix+'.scamp.xml'
;;       config.xsl_url = xslscamp
;;
;;       maxiter = 1
;;
;;       t0 = systime(1)
;;       for iter = 0L, maxiter do begin 
;;
;;          config.distort_degrees = '3'
;;          config.mosaic_type = 'UNCHANGED' ; 'LOOSE'
;;          config.pixscale_maxerr = '1.1'
;;          config.position_maxerr = '0.5'
;;          config.posangle_maxerr = '1.0'
;;          config.aheader_suffix = '.head'
;;
;;          im_scamp, catlist, config, configfile=configfile, silent=silent
;;
;;       endfor 
;;       splog, 'Total time to run scamp = ', (systime(1)-t0)/60.0, ' minutes.'
;;
;;    endif
       
return
end
