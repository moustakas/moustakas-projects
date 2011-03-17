pro flamingos_mosaics, sextractor=sextractor, scamp=scamp, swarp=swarp, test=test, jpeg=jpeg
; jm07aug20 - combine Anthony's "a" and "b" mosaics

    sexpath = sg1120_path(/sex)
    mosaicpath = flamingos_path(/mosaics)
    datapath = flamingos_path(/mosaics)

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
       imagelist = file_search(datapath+'sg1120[a,b]_ks.fits',count=nimage)
       weightlist = repstr(imagelist,'.fits','.weight.fits')

       catlist = datapath+file_basename(repstr(imagelist,'.fits','.cat'))
       radecreglist = repstr(catlist,'.cat','.radec.reg')
       seglist = repstr(catlist,'.cat','.seg.fits')

; initialize the SE configuration parameters
       config = init_sex_config(nimage)
       configfile = datapath+'sex.config'

       config.catalog_name = catlist
       config.weight_image = weightlist
;      config.flag_image = flaglist
       config.parameters_name = sexparam
       config.filter_name = sexconv
       config.starnnw_name = sexnnw

       config.catalog_type = 'FITS_LDAC'
       config.detect_thresh = 1.5
       config.analysis_thresh = 1.5
       config.weight_type = 'MAP_WEIGHT'
       config.weight_gain = 'N'
       config.interp_type = 'NONE'
       config.nthreads = 4

       config.seeing_fwhm = '0.8'

;      magzero = strtrim(sxpar(headfits(imagelist[ii]),'ZPT_FAT'),2)
       config.mag_zeropoint = '22.0'
       
       config.checkimage_type = 'NONE' ; SEGMENTATION
       config.checkimage_name = seglist

       mwrfits, config, configfile+'.fits', /create

; do it!
       t0 = systime(1)
       im_sex, imagelist, config, silent=silent ; do not pass CONFIGFILE
       splog, 'Total time to generate SE catalogs = ', $
         (systime(1)-t0)/60.0, ' minutes'

; make a QAplot and write out region files       
       sexqaplot = mosaicpath+'qaplot_Ks_sex.ps'
       dfpsplot, sexqaplot
       for ii = 0L, nimage-1L do begin
          cat = mrdfits(catlist[ii],2,/silent)
          xr = [8,22.5] & yrange1 = [-0.02,1.0] & yrange2 = [0.0,4.9]
          plot, [0], [0], /nodata, position=[0.14,0.55,0.96,0.93], xsty=1, ysty=1, $
            xrange=xr, yrange=yrange1, xtitle='', xtickname=replicate(' ',10), $
            ytitle='Class Star', charsize=1.8, xthick=5.0, $
            ythick=5.0, charthick=5.0, title=file_basename(repstr(imagelist[ii],'.fits',''))
          djs_oplot, cat.mag_auto, cat.class_star, psym=symcat(16)
          djs_plot, [0], [0], /nodata, /noerase, position=[0.14,0.12,0.96,0.55], $
            xsty=1, ysty=1, xrange=xr, yrange=yrange2, xtitle='Instrumental Magnitude', $
            ytitle='r_{eff} (arcsec)', charsize=1.8, xthick=5.0, $
            ythick=5.0, charthick=5.0
          djs_oplot, cat.mag_auto, cat.flux_radius*cat.awin_image*0.188, $
            psym=symcat(16), symsize=0.8

          splog, 'Witing DS9 region file '+file_basename(radecreglist[ii])
          write_ds9_regionfile, cat.xwin_world, cat.ywin_world, $
            filename=radecreglist[ii], color='red'
       endfor
       dfpsclose

;      t0 = systime(1)
;      for ii = 0L, n_elements(imagelist)-1L do begin
;         print, 'SExtracting '+imagelist[ii]
;         magzero = strtrim(sxpar(headfits(imagelist[ii]),'ZPT_FAT'),2)
;         spawn, 'sex '+imagelist[ii]+' -c '+sexconfig+' -CATALOG_NAME '+catlist[ii]+' -CATALOG_TYPE FITS_LDAC'+$
;           ' -PARAMETERS_NAME '+sexparam+' -DETECT_THRESH 1.5 -ANALYSIS_THRESH 1.5 -FILTER_NAME '+$
;           sexconv+' -STARNNW_NAME '+sexnnw+$
;           ' -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE '+weightlist[ii]+' -WEIGHT_THRESH 0 -SEEING_FWHM 0.5'+$
;           ' -MAG_ZEROPOINT '+magzero+' -NTHREADS 4 -CHECKIMAGE_TYPE NONE -CHECKIMAGE_NAME '+seglist[ii], /sh
;           ' -CHECKIMAGE_TYPE SEGMENTATION -CHECKIMAGE_NAME '+seglist[ii], /sh
;         cat = mrdfits(catlist[ii],2,/silent)
;         struct_print, struct_trimtags(cat,sele=['xwin_image','ywin_image'],$
;           format=['G15.8','G15.8']),file=xyreglist[ii], /no_head
;      endfor
;      splog, 'Total time to generate SE catalogs = ', $
;        (systime(1)-t0)/60.0, ' minutes.'
    endif 

; ---------------------------------------------------------------------------
; SCAMP

    if keyword_set(scamp) then begin
       catlist = file_search(datapath+'sg1120[a,b]_ks.cat')

; initialize the scamp configuration parameters
       config = init_scamp_config()
       configfile = mosaicpath+'scamp.config'

       config.astref_catalog = 'FILE'
       config.astrefcat_name = sexpath+'sg1120_sdss_dr7_refcat.cat'
;      config.astrefcat_name = sexpath+'sg1120_usnob_refcat.cat'
       config.astrefcent_keys = 'RA,DEC'
       config.astreferr_keys = 'ERR_A,ERR_B'
       config.astrefmag_key = 'MAG'

; boost the weights
;      config.astref_weight = '100.0'
;      config.crossid_radius = 1.5 ; this is important

       config.astrinstru_key = 'FILTER'
       config.photinstru_key = 'FILTER'
       config.magzero_key = 'ZPT_FAT'
       config.extinct_key = 'PHOT_K'
       config.photomflag_key = 'PHOTFLAG'
       config.save_refcatalog = 'N'
       config.refout_catpath = mosaicpath
       config.mergedoutcat_type = 'NONE'

       config.checkplot_type = strjoin(['ASTR_CHI2','ASTR_INTERROR1D','ASTR_INTERROR2D',$
         'ASTR_REFERROR1D','ASTR_REFERROR2D','DISTORTION','FGROUPS','PHOT_ERROR','PHOT_ZPCORR',$
         'PHOT_ZPCORR3D'],',')
       config.checkplot_name = strjoin('flamingos_'+['astr_chi2','astr_interror1d',$
         'astr_interror2d','astr_referror1d','referror2d','distort','fgroups',$
         'psphot_error','phot_zpcorr','phot_zpcorr3d'],',')
;      config.checkplot_type = 'NONE'
       config.checkplot_dev = 'JPEG'

       config.xml_name = mosaicpath+'flamingos.scamp.xml'
       config.xsl_url = xslscamp

       t0 = systime(1)
       maxiter = 1

; think carefully before increasing DEGREE      
       for iter = 0, maxiter do begin 
          case iter of
             0: begin
                config.distort_degrees = '1'
                config.mosaic_type = 'LOOSE'
                config.pixscale_maxerr = '1.1'
                config.position_maxerr = '5.0'
                config.posangle_maxerr = '3.0'
                config.aheader_suffix = '.ahead'
             end
             1: begin
                config.distort_degrees = '1,1'
                config.mosaic_type = 'LOOSE'
                config.pixscale_maxerr = '1.1'
                config.position_maxerr = '3.0'
                config.posangle_maxerr = '1.0'
                config.aheader_suffix = '.head'
             end
             else: begin
                config.distort_degrees = '1,1' ; '4,4'
                config.mosaic_type = 'LOOSE'
                config.pixscale_maxerr = '1.1'
                config.position_maxerr = '1.0'
                config.posangle_maxerr = '1.0'
                config.aheader_suffix = '.head'
             end
          endcase

          mwrfits, config, configfile+'.fits', /create
          im_scamp, catlist, config, configfile=configfile, silent=silent

       endfor 
       splog, 'Total time to run scamp = ', (systime(1)-t0)/60.0, ' minutes.'

       
;      astref_catalog = 'USNO-B1'
;      astref_band = 'Bj'
;      astrinstru_key = 'FILTER'
;      photinstru_key = 'FILTER'
;      magzero_key = 'ZPT_FAT'
;      extinct_key = 'PHOT_K'
;      photomflag_key = 'PHOTFLAG'
;      refout_catpath = mosaicpath
;
;      checkplot_type = ['DISTORTION','ASTR_CHI2','PHOT_ERROR','PHOT_ZPCORR3D','ASTR_INTERROR2D',$
;        'ASTR_INTERROR1D','ASTR_REFERROR2D','ASTR_REFERROR1D','PHOT_ZPCORR','ASTR_PIXERROR1D',$
;        'ASTR_REFSYSMAP']
;      checkplot_name = 'flamingos_'+['distort','astr_chi2','phot_error','phot_zpcorr3d',+$
;        'astr_interror2d','astr_interror1d','astr_referror2d','astr_referror1d',$
;        'phot_zpcorr','astr_pixerror1d','astr_refsysmap']
;      checkplot_type = strjoin(checkplot_type,',')
;      checkplot_name = strjoin(mosaicpath+checkplot_name,',')
;
;      t0 = systime(1)
;      for iter = 0L, 2L do begin
;
;         case iter of
;            0L: begin
;               degree = '1'
;               mosaic_type = 'LOOSE'
;               position_maxerr = '1.0'
;               posangle_maxerr = '1.0'
;               aheader_suffix = '.ahead'
;            end
;            1L: begin
;               degree = '1'
;               mosaic_type = 'LOOSE'
;               position_maxerr = '0.5'
;               posangle_maxerr = '1.0'
;               aheader_suffix = '.head'
;            end
;            else: begin
;               degree = '1'
;               mosaic_type = 'LOOSE'
;               position_maxerr = '0.5'
;               posangle_maxerr = '1.0'
;               aheader_suffix = '.head'
;            end
;         endcase
;
;         spawn, 'scamp '+strjoin(catlist,' ')+' -c '+scampconfig+' -MERGEDOUTCAT_TYPE NONE'+$ ;FITS_LDAC'+$
;           ' -PIXSCALE_MAXERR 1.1 -POSANGLE_MAXERR '+posangle_maxerr+' -POSITION_MAXERR '+position_maxerr+$
;           ' -ASTREF_CATALOG '+astref_catalog+' -ASTREF_BAND '+astref_band+' -SAVE_REFCATALOG Y -REFOUT_CATPATH '+refout_catpath+$
;           ' -CHECKPLOT_DEV PSC -CHECKPLOT_TYPE '+checkplot_type+' -CHECKPLOT_NAME '+checkplot_name+$
;           ' -AHEADER_SUFFIX '+aheader_suffix+' -DISTORT_DEGREES '+degree+' -MOSAIC_TYPE '+mosaic_type+$
;           ' -ASTRINSTRU_KEY '+astrinstru_key+' -SOLVE_PHOTOM Y -PHOTINSTRU_KEY '+photinstru_key+$
;           ' -MAGZERO_KEY '+magzero_key+' -EXTINCT_KEY '+extinct_key+' -PHOTOMFLAG_KEY '+photomflag_key+$
;           ' -NTHREADS 4 -WRITE_XML N -VERBOSE_TYPE NORMAL', /sh
;
;      endfor 
;      splog, 'Total time to run scamp = ', (systime(1)-t0)/60.0, ' minutes.'

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
          config.pixel_scale = '0.317'
          config.center = '11:20:22.0,-12:04:08.0'
          config.image_size = '3200,3890'
       endelse 
       
       config.blank_badpixels = 'Y'
       config.interpolate = 'N'
       config.write_fileinfo = 'Y'
       config.celestial_type = 'EQUATORIAL'
       config.copy_keywords = 'OBJECT,FILTER'
       config.xml_name = mosaicpath+'flamingos.swarp.xml'
       config.xsl_url = xslswarp

; build the Ks-band mosaic
       
       config.write_xml = 'Y'
       config.subtract_back = 'Y'
       config.resampling_type = 'LANCZOS3'
       config.weight_type = 'MAP_WEIGHT'

       config.combine_type = 'WEIGHTED'
       
       configfile = mosaicpath+'swarp.config'

       config.imageout_name = mosaicpath+'sg1120_Ks.fits'
       config.weightout_name = repstr(config.imageout_name,'.fits','.weight.fits')

       imagelist = file_search(datapath+'sg1120[a,b]_ks.fits')
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
               sxaddpar, hdr, 'OBJECT', 'SG1120 Ks'
               modfits, fixme[ff], 0, hdr
           endfor
       endif
          
   endif           

stop   
   
; ---------------------------------------------------------------------------
; generate a greyscale TIFF and JPEG
    if keyword_set(jpeg) then begin

       fitsfile = mosaicpath+'sg1120_ks.fits'
       tifffile = mosaicpath+'sg1120_ks.tiff'
       jpegfile = mosaicpath+'sg1120_ks.jpeg'

       splog, 'Writing '+tifffile
       spawn, 'stiff -c '+stiffconfig+' '+fitsfile+' -OUTFILE_NAME '+tifffile+$
         ' -BINNING 2 -GAMMA_FAC 1.0'

       splog, 'Writing '+jpegfile
       spawn, '/usr/bin/convert '+tifffile+' '+jpegfile
       
    endif

return
end
