pro hst_mosaics, sextractor=sextractor, scamp=scamp, swarp=swarp, $
  small=small, test=test, blakeslee=blakeslee, jpeg=jpeg
; jm07aug - generate mosaics from the individual HST pointings;
;           running SEXTRACTOR and SCAMP doesn't really work because
;           there are too few USNO stars; instead, build the mosaic
;           using SWARP, and then register the full mosaic using
;           SG1120_REGISTER_MOSAICS 

    sexpath = sg1120_path(/sex)
    mosaicpath = hst_path(/mosaics)
    datapath = hst_path()
    if keyword_set(blakeslee) then begin
       datapath = datapath+'blakeslee_drizzle/'
       weightsuffix = '_weight'
    endif else begin
       datapath = datapath+'gonzalez_drizzle/'
       weightsuffix = '.weight'
    endelse

    stiffconfig = sexpath+'default.stiff'
    scampconfig = sexpath+'default.scamp'
    swarpconfig = sexpath+'default.swarp'
    sexconfig = sexpath+'default.sex'
    sexparam = sexpath+'sg1120.sex.param'
    sexconv = sexpath+'default.conv'
    sexnnw = sexpath+'default.nnw'

; ---------------------------------------------------------------------------    
; generate SE catalogs
; ---------------------------------------------------------------------------    
    
    if keyword_set(sextractor) then begin

       splog, 'Are you sure you want to do this?' & stop
       
       if keyword_set(blakeslee) then $
         imagelist = file_search(datapath+'supergroup-pos*_f814w_drz_sci.fits') else $
         imagelist = file_search(datapath+'hst_f814w_drz_pos??.fits')

       weightlist = repstr(imagelist,'.fits',weightsuffix+'.fits')
       catlist = datapath+file_basename(repstr(imagelist,'.fits','.cat'))
;      catlist = mosaicpath+file_basename(repstr(imagelist,'.fits','.cat'))
       seglist = repstr(catlist,'.cat','.seg.fits')
       xyreglist = repstr(catlist,'.cat','.xy.reg')
       radecreglist = repstr(catlist,'.cat','.radec.reg')

       magzero = string(25.501+2.5*alog10(2000.0),format='(F7.4)')

       t0 = systime(1)
       for ii = 0L, n_elements(imagelist)-1L do begin
          print, 'SExtracting '+imagelist[ii]
          spawn, 'sex '+imagelist[ii]+' -c '+sexconfig+' -CATALOG_NAME '+catlist[ii]+' -CATALOG_TYPE FITS_LDAC'+$
            ' -PARAMETERS_NAME '+sexparam+' -DETECT_THRESH 2.0 -ANALYSIS_THRESH 2.0 -FILTER_NAME '+sexconv+' -STARNNW_NAME '+sexnnw+$
            ' -DETECT_MINAREA 25 -DEBLEND_NTHRESH 16 -DEBLEND_MINCONT 0.015'+$
            ' -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE '+weightlist[ii]+' -WEIGHT_THRESH 0,1E29 -SEEING_FWHM 0.105'+$
            ' -MAG_ZEROPOINT '+magzero+' -NTHREADS 4 -CHECKIMAGE_TYPE NONE -CHECKIMAGE_NAME '+seglist[ii], /sh
;           ' -CHECKIMAGE_TYPE SEGMENTATION -CHECKIMAGE_NAME '+seglist[ii], /sh
          cat = mrdfits(catlist[ii],2,/silent)
;         struct_print, struct_trimtags(cat,sele=['xwin_image','ywin_image'],$
;           format=['G15.8','G15.8']),file=xyreglist[ii], /no_head
;         struct_print, struct_trimtags(cat,sele=['xwin_world','ywin_world'],$
;           format=['F15.8','F15.8']),file=radecreglist[ii], /no_head
       endfor
       splog, 'Total time to generate SE catalogs = ', (systime(1)-t0)/60.0, ' minutes.'

    endif

; ---------------------------------------------------------------------------
; SCAMP
; ---------------------------------------------------------------------------

    if keyword_set(scamp) then begin

       splog, 'Are you sure you want to do this?' & stop

       catlist = file_search(mosaicpath+'hst_f814w_drz_pos??.cat')

       astref_catalog = 'USNO-B1'
       astref_band = 'Bj'
       astrinstru_key = 'FILTER'
       photinstru_key = 'FILTER'
       magzero_key = 'MAGZERO'
       extinct_key = 'PHOT_K'
       photomflag_key = 'PHOTFLAG'
       refout_catpath = mosaicpath

       checkplot_type = ['DISTORTION','ASTR_CHI2','PHOT_ERROR','PHOT_ZPCORR3D','ASTR_INTERROR2D',$
         'ASTR_INTERROR1D','ASTR_REFERROR2D','ASTR_REFERROR1D','PHOT_ZPCORR','ASTR_PIXERROR1D',$
         'ASTR_REFSYSMAP']
       checkplot_name = 'vimos'+suffix+'_'+['distort','astr_chi2','phot_error','phot_zpcorr3d',+$
         'astr_interror2d','astr_interror1d','astr_referror2d','astr_referror1d',$
         'phot_zpcorr','astr_pixerror1d','astr_refsysmap']
       checkplot_type = strjoin(checkplot_type,',')
       checkplot_name = strjoin(mosaicpath+checkplot_name,',')

       t0 = systime(1)
       for iter = 0L, 5L do begin ; three iterations is good

          case iter of
             0L: begin
                degree = '1'
                mosaic_type = 'LOOSE'
                position_maxerr = '1.0'
                posangle_maxerr = '1.0'
                aheader_suffix = '.ahead'
             end
             1L: begin
                degree = '1'
                mosaic_type = 'LOOSE'
;               mosaic_type = 'UNCHANGED'
                position_maxerr = '1.0'
                posangle_maxerr = '1.0'
                aheader_suffix = '.head'
             end
             else: begin
                degree = '1'
                mosaic_type = 'LOOSE'
;               mosaic_type = 'UNCHANGED'
                position_maxerr = '1.0'
                posangle_maxerr = '1.0'
                aheader_suffix = '.head'
             end
          endcase

          spawn, 'scamp '+strjoin(catlist,' ')+' -c '+scampconfig+' -MERGEDOUTCAT_TYPE NONE'+$ ;FITS_LDAC'+$
            ' -PIXSCALE_MAXERR 1.1 -POSANGLE_MAXERR '+posangle_maxerr+' -POSITION_MAXERR '+position_maxerr+$
            ' -ASTREF_CATALOG '+astref_catalog+' -ASTREF_BAND '+astref_band+' -SAVE_REFCATALOG Y -REFOUT_CATPATH '+refout_catpath+$
            ' -CHECKPLOT_DEV PSC -CHECKPLOT_TYPE '+checkplot_type+' -CHECKPLOT_NAME '+checkplot_name+$
            ' -AHEADER_SUFFIX '+aheader_suffix+' -DISTORT_DEGREES '+degree+' -MOSAIC_TYPE '+mosaic_type+$
;           ' -CROSSID_RADIUS 3.0 -FWHM_THRESHOLDS 0.0,100.0'+$; -STABILITY_TYPE EXPOSURE'+$
            ' -ASTRINSTRU_KEY '+astrinstru_key+' -NTHREADS 1 -WRITE_XML N -VERBOSE_TYPE NORMAL', /sh

       endfor 
       splog, 'Total time to run scamp = ', (systime(1)-t0)/60.0, ' minutes.'

    endif
    
; ---------------------------------------------------------------------------
; SWARP
; ---------------------------------------------------------------------------

    if keyword_set(swarp) then begin
       
       if keyword_set(test) then begin
          radec = '00:00:00.0,+00:00:00.0'
          center_type = 'ALL' 
          pixelscale_type = 'MANUAL' 
          if keyword_set(blakeslee) then begin
             if keyword_set(small) then pixel_scale = '0.10' else pixel_scale = '0.035'
             suffix = '_blakeslee_test_'+pixel_scale
          endif else begin
             if keyword_set(small) then pixel_scale = '0.10' else pixel_scale = '0.03'
             suffix = '_test_'+pixel_scale
          endelse
          image_size = '0'
          header_only = 'Y'
       endif else begin
          radec = '11:20:16.0,-12:03:36.0'
          center_type = 'MANUAL' 
          pixelscale_type = 'MANUAL' 
          if keyword_set(blakeslee) then begin
             if keyword_set(small) then begin
                pixel_scale = '0.10' 
                image_size = '7470,11600'
             endif else begin
                pixel_scale = '0.035' 
                image_size = '21380,33150'
             endelse
             suffix = '_blakeslee_'+pixel_scale
          endif else begin
             if keyword_set(small) then begin
                pixel_scale = '0.10' 
                image_size = '7470,11600'
             endif else begin
                pixel_scale = '0.03'
                image_size = '24900,38660'
             endelse 
             suffix = '_'+pixel_scale
          endelse
          header_only = 'N'
       endelse

       combinetype = 'WEIGHTED'
       gain_keyword = 'JUNK'
;      gain_keyword = 'CCDGAIN'
       projection_err = '0.001'
       keywords = 'JUNK'
;      keywords = 'OBJECT,FILTER'
       
       mosaic_file = mosaicpath+'sg1120_hst'+suffix+'.fits'
       mosaic_weightfile = repstr(mosaic_file,'.fits','.weight.fits')
       
; for some reason if POS 10 appears first then I get an error in
; swarp, so order the files as POS 1-10
       
       if keyword_set(blakeslee) then $
         imagelist = [file_search(datapath+'supergroup-pos?_f814w_drz_sci.fits'),$
         file_search(datapath+'supergroup-pos10_f814w_drz_sci.fits')] else $
         imagelist = [file_search(datapath+'hst_f814w_drz_pos0?.fits'),file_search(datapath+'hst_f814w_drz_pos10.fits')]

       splog, 'Building '+mosaic_file+' from '+string(n_elements(imagelist),format='(I0)')+' images.'
       spawn, 'swarp '+strjoin(imagelist,',')+' -c '+swarpconfig+' -HEADER_ONLY '+header_only+$
         ' -IMAGEOUT_NAME '+mosaic_file+' -WEIGHTOUT_NAME '+mosaic_weightfile+$
         ' -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_SUFFIX '+weightsuffix+'.fits -WEIGHT_THRESH 0'+$
         ' -COMBINE_TYPE '+combinetype+' -BLANK_BADPIXELS Y -CELESTIAL_TYPE EQUATORIAL'+$
         ' -PROJECTION_TYPE TAN -PROJECTION_ERR '+projection_err+' -CENTER_TYPE '+center_type+' -CENTER '+radec+$
         ' -PIXELSCALE_TYPE '+pixelscale_type+' -PIXEL_SCALE '+pixel_scale+' -IMAGE_SIZE '+image_size+$
         ' -RESAMPLING_TYPE LANCZOS3 -INTERPOLATE N -GAIN_KEYWORD '+gain_keyword+' -SUBTRACT_BACK N'+$
         ' -COPY_KEYWORDS '+keywords+$
         ' -COMBINE_BUFSIZE 1024 -WRITE_FILEINFO N -VERBOSE_TYPE NORMAL -NTHREADS 4', /sh

       if (not keyword_set(test)) then begin
          hdr = headfits(mosaic_file)
          sxaddpar, hdr, 'OBJECT', 'SG1120'
          sxaddpar, hdr, 'FILTER', 'F814W'
          modfits, mosaic_file, 0, hdr

          hdr = headfits(mosaic_weightfile)
          sxaddpar, hdr, 'OBJECT', 'SG1120'
          sxaddpar, hdr, 'FILTER', 'F814W'
          modfits, mosaic_weightfile, 0, hdr
       endif
          
    endif
       
; ---------------------------------------------------------------------------
; generate a greyscale TIFF and JPEG
; ---------------------------------------------------------------------------

    if keyword_set(jpeg) then begin

       if keyword_set(blakeslee) then suffix = '_blakeslee_0.10' else suffix = '_0.10'
       
       fitsfile = mosaicpath+'sg1120_hst'+suffix+'.fits'
       tifffile = mosaicpath+'sg1120_hst'+suffix+'.tiff'
       jpegfile = mosaicpath+'sg1120_hst'+suffix+'.jpeg'

       splog, 'Writing '+tifffile
       spawn, 'stiff -c '+stiffconfig+' '+fitsfile+' -OUTFILE_NAME '+tifffile+$
         ' -BINNING 2'

       splog, 'Writing '+jpegfile
       spawn, '/usr/bin/convert '+tifffile+' '+jpegfile
       
    endif

return
end
