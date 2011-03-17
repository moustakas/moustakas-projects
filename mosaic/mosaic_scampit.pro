pro mosaic_scampit, preproc=preproc, badpix=badpix, sextractor=sextractor, scamp=scamp, swarp=swarp, test=test, jpeg=jpeg
; jm07jan25nyu - generate all the mosaics

; rsync -auv line:"/global/data/scr/morad/4meter/f/fobj12[5-9].fits" /Volumes/WDexternal/data/mosaic/rawdata/
; rsync -auv line:"/global/data/scr/morad/4meter/f/fobj130.fits"     /Volumes/WDexternal/data/mosaic/rawdata/
; rsync -auv line:"/global/data/scr/morad/4meter/f/fobj13[2-6].fits" /Volumes/WDexternal/data/mosaic/rawdata/
; rsync -auv line:"/global/data/scr/morad/4meter/f/fobj13[8-9].fits" /Volumes/WDexternal/data/mosaic/rawdata/
; rsync -auv line:"/global/data/scr/morad/4meter/f/fobj14[0-2].fits" /Volumes/WDexternal/data/mosaic/rawdata/
; rsync -auv line:"/global/data/scr/morad/4meter/f/fobj179.fits"     /Volumes/WDexternal/data/mosaic/rawdata/
; rsync -auv line:"/global/data/scr/morad/4meter/f/fobj18[0-3].fits" /Volumes/WDexternal/data/mosaic/rawdata/
; rsync -auv line:"/global/data/scr/morad/4meter/f/fobj22[6-9].fits" /Volumes/WDexternal/data/mosaic/rawdata/
; rsync -auv line:"/global/data/scr/morad/4meter/f/fobj23[0-2].fits" /Volumes/WDexternal/data/mosaic/rawdata/
    
;   path = '/Volumes/WDexternal/data/mosaic/'
    path = '/global/bias1/ioannis/mosaic/'
    rootpath = '/global/data/scr/morad/4meter/f/'
    datapath = path+'data_chip7/'
;   datapath = path+'data/'
    sexpath = path+'sex/'
    mosaicpath = path+'mosaics/'
    
    scampconfig = sexpath+'mosaic.scamp'
    swarpconfig = sexpath+'mosaic.swarp'
    sexconfig = sexpath+'mosaic.sex'
    sexparam = sexpath+'mosaic.sex.param'
    sexconv = sexpath+'default.conv'
    sexnnw = sexpath+'default.nnw'

; ---------------------------------------------------------------------------    

    if keyword_set(preproc) then begin
;      imagelist = file_search(rootpath+'fobj???.fits')
;      imagelist = file_search(rootpath+'fobj12[5-7].fits')
       imagelist = datapath+'fobj'+['125','126','127','128','129','130','138','139','140','141','142',$
         '132','133','134','135','136','179','180','181','182','183',$
         '226','227','228','229','230','231','232']+'.fits'
       rawimagelist = rootpath+file_basename(imagelist)
       newimagelist = datapath+file_basename(imagelist)
       t0 = systime(1)
       for ii = 0L, n_elements(imagelist)-1L do begin ; strip the TNX header entries
          bighdr = headfits(rawimagelist[ii],ext=0)
          pixscale = sxpar(bighdr,'PIXSCAL1') ; [arcsec/pixel]
          mwrfits, 0, newimagelist[ii], /create
;         for iext = 7L, 7L do begin
          for iext = 1L, 8L do begin
             splog, 'Processing '+imagelist[ii]+', extension '+string(iext,format='(I0)')
             im = mrdfits(rawimagelist[ii],iext,hdr,/silent)
             astr = hogg_make_astr(sxpar(hdr,'CRVAL1'),sxpar(hdr,'CRVAL2'),$
               sxpar(hdr,'NAXIS1')*pixscale/3600.0D,sxpar(hdr,'NAXIS2')*pixscale/3600.0D,$
               orientation=90.0,pixscale=pixscale/3600.0D)
             astr.crpix = [sxpar(hdr,'CRPIX1'),sxpar(hdr,'CRPIX2')] ; NOTE!
             mkhdr, newhdr, im  ; generate a basic FITS header
             sxdelpar, newhdr, 'COMMENT'
             sxdelpar, newhdr, 'DATE'
             sxaddpar, newhdr, 'DATE-OBS', sxpar(hdr,'DATE-OBS')
             sxaddpar, newhdr, 'OBJECT', sxpar(hdr,'OBJECT')
             sxaddpar, newhdr, 'GAIN', float(sxpar(hdr,'GAIN'))
             sxaddpar, newhdr, 'RDNOISE', float(sxpar(hdr,'RDNOISE'))
             sxaddpar, newhdr, 'EXPTIME', float(sxpar(hdr,'EXPTIME'))
             sxaddpar, newhdr, 'AIRMASS', float(sxpar(hdr,'AIRMASS'))
             sxaddpar, newhdr, 'FILTER', sxpar(hdr,'FILTER')
             sxdelpar, newhdr, 'HISTORY'
             putast, newhdr, astr
             splog, 'Writing '+newimagelist[ii]+', extension '+string(iext,format='(I0)')
             mwrfits, im, newimagelist[ii], newhdr
          endfor
          print
       endfor
       splog, 'Total time to pre-process images = ', (systime(1)-t0)/60.0, ' minutes'
       
    endif

; ---------------------------------------------------------------------------    

    if keyword_set(badpix) then begin
       masklist = path+'mosaic_bitmask_'+['g','r','i']+'.fits'
       objectlist = ['g','r','i']+'-band badpix mask'
       outmasklist = datapath+file_basename(repstr(masklist,'bitmask','badpix'))
       for ii = 0L, 2L do begin
          mwrfits, 0, outmasklist[ii], headfits(masklist[ii]), /create
          for iext = 7L, 7L do begin
;         for iext = 1L, 8L do begin
             im = mrdfits(masklist[ii],iext,hdr)
             sxaddpar, hdr, 'OBJECT', objectlist[ii]
             newim = im eq 0
             mwrfits, newim, outmasklist[ii], hdr
          endfor
       endfor
    endif

; ---------------------------------------------------------------------------    

    if keyword_set(sextractor) then begin

;      imagelist = datapath+'fobj'+['125','126','127']+'.fits'
;      imagelist = datapath+'fobj'+['125','126','127','128','129','130','138','139','140','141','142']+'.fits'
       imagelist = file_search(datapath+'fobj???.fits')
       weightlist = repstr(imagelist,'.fits','.weight.fits')
       catlist = repstr(imagelist,'.fits','.cat')
       skylist = repstr(imagelist,'.fits','.sky.fits')

       t0 = systime(1)
;      for ii = 0L, 0L do begin
       for ii = 0L, n_elements(imagelist)-1L do begin

          if strmatch(sxpar(headfits(imagelist[ii],ext=1),'FILTER'),'*g SDSS*') then weight_image = datapath+'mosaic_badpix_g.fits'
          if strmatch(sxpar(headfits(imagelist[ii],ext=1),'FILTER'),'*r SDSS*') then weight_image = datapath+'mosaic_badpix_r.fits'
          if strmatch(sxpar(headfits(imagelist[ii],ext=1),'FILTER'),'*i SDSS*') then weight_image = datapath+'mosaic_badpix_i.fits'

          print, 'SExtracting '+imagelist[ii]
          spawn, 'sex '+imagelist[ii]+' -c '+sexconfig+' -CATALOG_NAME '+catlist[ii]+' -CATALOG_TYPE FITS_LDAC'+$
            ' -PARAMETERS_NAME '+sexparam+' -DETECT_THRESH 1.0 -ANALYSIS_THRESH 1.0 -FILTER_NAME '+sexconv+' -STARNNW_NAME '+sexnnw+$
            ' -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE '+weight_image+' -WEIGHT_THRESH 0'+$
;           ' -CHECKIMAGE_TYPE BACKGROUND -CHECKIMAGE_NAME '+skylist[ii]+$
            ' -MEMORY_BUFSIZE 8192 -VERBOSE_TYPE NORMAL -NTHREADS 4', /sh

       endfor
       splog, 'Total time to run = ', (systime(1)-t0)/60.0, ' minutes.'

    endif

; ---------------------------------------------------------------------------    

    if keyword_set(scamp) then begin

       catlist = file_search(datapath+'fobj???.cat')
;      catlist = datapath+'fobj'+['125','126','127','128','129','130','138','139','140','141','142']+'.cat'
;      catlist = datapath+'fobj'+['125','126','127','128','129','130','138','139','140','141','142',$
;                                 '132','133','134','135','136','179','180','181','182','183',$
;                                 '226','227','228','229','230','231','232']+'.cat'
;      catlist = file_search(datapath+'fobj12[5-7].cat')

;      astref_catalog = 'USNO-B1'
;      astref_catalog = 'SDSS-R3'
       astref_catalog = 'SDSS-R5'
       astrinstru_key = 'FILTER'
       photinstru_key = 'FILTER'
       magzero_key = 'PHOT_C'
       extinct_key = 'PHOT_K'
       photomflag_key = 'PHOTFLAG'

;      for iter = 0L, 3L do begin
       for iter = 0L, 2L do begin

          case iter of
             0L: begin
                degree = '3'
                mosaic_type = 'UNCHANGED'
                crossid_radius = '7.0'
                position_maxerr = '2.0'
                posangle_maxerr = '3.0'
                aheader_suffix = '.ahead'
             end
              1L: begin
                degree = '3'
;               mosaic_type = 'UNCHANGED'
;               mosaic_type = 'LOOSE'
                mosaic_type = 'FIX_FOCALPLANE'
                crossid_radius = '5.0'
                position_maxerr = '0.5'
                posangle_maxerr = '2.0'
                aheader_suffix = '.head'
             end
             else: begin
                degree = '3'
;               mosaic_type = 'UNCHANGED'
;               mosaic_type = 'LOOSE'
                mosaic_type = 'FIX_FOCALPLANE'
                crossid_radius = '3.0'
                position_maxerr = '0.1'
                posangle_maxerr = '1.0'
                aheader_suffix = '.head'
             end
          endcase

          spawn, 'scamp '+strjoin(catlist,' ')+' -c '+scampconfig+' -CHECKPLOT_DEV PSC -NTHREADS 4'+$
            ' -PIXSCALE_MAXERR 1.1 -POSANGLE_MAXERR '+posangle_maxerr+' -POSITION_MAXERR '+position_maxerr+$
            ' -ASTREF_CATALOG '+astref_catalog+' -MERGEDOUTCAT_TYPE FITS_LDAC -SAVE_REFCATALOG Y'+$
            ' -AHEADER_SUFFIX '+aheader_suffix+' -DISTORT_DEGREES '+degree+' -MOSAIC_TYPE '+mosaic_type+$
            ' -CROSSID_RADIUS '+crossid_radius+$; -STABILITY_TYPE EXPOSURE'+$;'-SN_THRESHOLDS 3.0,100.0'+$
            ' -ASTRINSTRU_KEY '+astrinstru_key+' -WRITE_XML N -VERBOSE_TYPE NORMAL', /sh

       endfor    
    
    endif

; ---------------------------------------------------------------------------    
    
    if keyword_set(swarp) then begin

; global parameters       

       combinetype = 'WEIGHTED'
       projection_err = '0.001'
       keywords = 'OBJECT,FILTER'
       pixelscale_type = 'MANUAL' 

       pixel_scale = '0.258'
       gain_keyword = 'GAIN'
       
; generate six mosaics (groups) in three different bandpasses

       bandpass = ['g','r','i']
       ngroup = 6L
       
       weight_image = datapath+'mosaic_badpix_'+bandpass+'.fits'
       
       for ig = 4L, 4L do begin
;      for ig = 1L, ngroup do begin
      
;         for ib = 0L, 0L do begin
          for ib = 0L, n_elements(bandpass)-1L do begin
          
             case ig of
                1L: begin       ; Group 1 - UMa
                   case ib of
                      0L: flist = ['059','060','061','062','063','185','186','187','188','189',$ ; scamp problem with 251
                        '234','235','236','237','238','239','240']
                      1L: flist = ['064','065','066','068','076','077','078','079','080',$
                        '191','192','193','194','242','243','244','245','246']
                      2L: flist = ['070','071','072','073','074','223','224','225',$ ; scamp problem with 251
                        '247','248','249','250','251']
                   endcase
                   object = 'UMa'
                   radec = '10:36:40.0,+51:57:40.0'
                   image_size = '17343,10253'
                end
                2L: begin
                   case ib of   ; Group 2 - Willman Candidate 1
                      0L: flist = ['082']
                      1L: flist = ['083']
                      2L: flist = '' ; no data
                   endcase
                   object = 'Willman_candidate_1'
                   radec = '13:54:46.7,+40:04:28.6'
                   image_size = '8718,8756'
                end
                3L: begin
                   case ib of   ; Group 3 - Willman Candidate 2
                      0L: flist = ['087','088']
                      1L: flist = ['089','090']
                      2L: flist = '' ; no data
                   endcase
                   object = 'Willman_candidate_2'
                   radec = '14:30:02.9,+48:11:38.4'
                   image_size = '8722,8758'
                end                
                4L: begin
                   case ib of   ; Group 4 - Willman 1
;                     0L: flist = ['125']
;                     0L: flist = ['125','126']
;                     0L: flist = ['125','126','127']
                      0L: flist = ['125','126','127','128','129','130','138','139','140','141','142']
                      1L: flist = ['132','133','134','135','136','179','180','181','182','183']
                      2L: flist = ['226','227','228','229','230','231','232']
                   endcase
                   object = 'Willman_1'
                   radec = '10:48:50.0,+51:02:30.0'
                   image_size = '5000,3000'
;                  radec = '10:49:45.0,+51:05:40.0'
;                  image_size = '9200,9200'
                end
                5L: begin
                   case ib of   ; Group 5 - Willman Candidate 3
                      0L: flist = ['145','146']
                      1L: flist = ['147','148']
                      2L: flist = '' ; no data
                   endcase
                   object = 'Willman_candidate_3'
                   radec = '15:28:29.0,+53:03:36.0'
                   image_size = '8701,8763'
                end
                6L: begin
                   case ib of   ; Group 6 - M51
                      0L: flist = ['149','152','155']
                      1L: flist = ['150','153','156']
                      2L: flist = ['151','154','157']
                   endcase
                   object = 'M51'
                   radec = '13:29:54.0,+47:11:50.0'
                   image_size = '11948,10462'
                end
             endcase
             
             catheadlist = file_search(datapath+'fobj'+flist+'.head',count=nobj)
             if (nobj gt 0L) then begin
             
                imagelist = repstr(catheadlist,'.head','.fits')
                mosaic_file = mosaicpath+object+'_'+bandpass[ib]+'.fits'
                mosaic_weightfile = repstr(mosaic_file,'.fits','.weight.fits')

; if testing, then over-write these parameters         

                if keyword_set(test) then begin
                   center_type = 'ALL' 
                   header_only = 'Y'
                   radec = '00:00:00.0,+00:00:00.0'
                   image_size = '0'
                endif else begin
                   center_type = 'MANUAL' 
                   header_only = 'N'
                   radec = radec
                   image_size = image_size
                endelse
                
                splog, 'Building '+mosaic_file+' from '+string(n_elements(imagelist),format='(I0)')+' images.'
                spawn, 'swarp '+strjoin(imagelist,',')+' -c '+swarpconfig+' -HEADER_ONLY '+header_only+$
                  ' -IMAGEOUT_NAME '+mosaic_file+' -WEIGHTOUT_NAME '+mosaic_weightfile+$
                  ' -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE '+weight_image[ib]+' -WEIGHT_THRESH 0'+$
                  ' -COMBINE_TYPE '+combinetype+' -BLANK_BADPIXELS Y -CELESTIAL_TYPE EQUATORIAL'+$
                  ' -PROJECTION_TYPE TAN -PROJECTION_ERR '+projection_err+' -CENTER_TYPE '+center_type+' -CENTER '+radec+$
                  ' -PIXELSCALE_TYPE '+pixelscale_type+' -PIXEL_SCALE '+pixel_scale+' -IMAGE_SIZE '+image_size+$
                  ' -RESAMPLING_TYPE LANCZOS3 -INTERPOLATE N -GAIN_KEYWORD '+gain_keyword+' -SUBTRACT_BACK Y'+$
                  ' -COMBINE_BUFSIZE 512 -COPY_KEYWORDS '+keywords+' -WRITE_FILEINFO N -VERBOSE_TYPE NORMAL -NTHREADS 4', /sh

             endif

          endfor ; close bandpass loop

       endfor ; close group loop

    endif

    if keyword_set(jpeg) then begin

; Willman 1
        jpegfile = mosaicpath+'Willman_1.jpeg'
        splog, 'Reading '+mosaicpath+'Willman_1_g.fits'
        gim = mrdfits(mosaicpath+'Willman_1_g.fits',0,ghdr,/silent)
        imsize = size(gim,/dim)
;       nx = 512L
        nx = 512 * (min(imsize)/512)
        gim_nw = (temporary(gim))$
          [imsize[0]/2L-nx/2L:imsize[0]/2L+nx/2L-1L, $
           imsize[1]/2L-nx/2L:imsize[1]/2L+nx/2L-1L]
        splog, 'Reading '+mosaicpath+'Willman_1_r.fits'
        rim_nw = (mrdfits(mosaicpath+'Willman_1_r.fits',0,rhdr,/silent))$
          [imsize[0]/2L-nx/2L:imsize[0]/2L+nx/2L-1L, $
           imsize[1]/2L-nx/2L:imsize[1]/2L+nx/2L-1L]
        splog, 'Reading '+mosaicpath+'Willman_1_i.fits'
        iim_nw = (mrdfits(mosaicpath+'Willman_1_i.fits',0,ihdr,/silent))$
          [imsize[0]/2L-nx/2L:imsize[0]/2L+nx/2L-1L, $
           imsize[1]/2L-nx/2L:imsize[1]/2L+nx/2L-1L]
        splog, 'Writing '+jpegfile
        scales= [0.3,1.0,1.5]*1000.0
        nw_rgb_make, iim_nw, rim_nw, gim_nw, name=jpegfile, $
          scales=scales, nonlinearity=3.0, rebinfactor=2, quality=95

stop        
        
; UMa
        jpegfile = mosaicpath+'UMa_gri'
        splog, 'Reading '+mosaicpath+'UMa_g.fits'
        gim = mrdfits(mosaicpath+'UMa_g.fits',0,ghdr,/silent)
        imsize = size(gim,/dim)
        nx = 512 * (min(imsize)/512)
        gim_nw = (temporary(gim))$
          [imsize[0]/2L-nx/2L:imsize[0]/2L+nx/2L-1L, $
           imsize[1]/2L-nx/2L:imsize[1]/2L+nx/2L-1L]
        splog, 'Reading '+mosaicpath+'UMa_r.fits'
        rim_nw = (mrdfits(mosaicpath+'UMa_r.fits',0,rhdr,/silent))$
          [imsize[0]/2L-nx/2L:imsize[0]/2L+nx/2L-1L, $
           imsize[1]/2L-nx/2L:imsize[1]/2L+nx/2L-1L]
        splog, 'Reading '+mosaicpath+'UMa_i.fits'
        iim_nw = (mrdfits(mosaicpath+'UMa_i.fits',0,ihdr,/silent))$
          [imsize[0]/2L-nx/2L:imsize[0]/2L+nx/2L-1L, $
           imsize[1]/2L-nx/2L:imsize[1]/2L+nx/2L-1L]
        splog, 'Writing '+jpegfile+'.jpg'
        nw_rgb_make, iim_nw, rim_nw, gim_nw, name=jpegfile, $
          scales=scales, nonlinearity=3.0, rebinfactor=2, quality=95
       
; M51       
        jpegfile = mosaicpath+'M51_gri'
        splog, 'Reading '+mosaicpath+'M51_g.fits'
        gim = mrdfits(mosaicpath+'M51_g.fits',0,ghdr,/silent)
        imsize = size(gim,/dim)
        nx = 512 * (min(imsize)/2048)
        gim_nw = (temporary(gim))$
          [imsize[0]/2L-nx/2L:imsize[0]/2L+nx/2L-1L, $
           imsize[1]/2L-nx/2L:imsize[1]/2L+nx/2L-1L]
        splog, 'Reading '+mosaicpath+'M51_r.fits'
        rim_nw = (mrdfits(mosaicpath+'M51_r.fits',0,rhdr,/silent))$
          [imsize[0]/2L-nx/2L:imsize[0]/2L+nx/2L-1L, $
           imsize[1]/2L-nx/2L:imsize[1]/2L+nx/2L-1L]
        splog, 'Reading '+mosaicpath+'M51_i.fits'
        iim_nw = (mrdfits(mosaicpath+'M51_i.fits',0,ihdr,/silent))$
          [imsize[0]/2L-nx/2L:imsize[0]/2L+nx/2L-1L, $
           imsize[1]/2L-nx/2L:imsize[1]/2L+nx/2L-1L]
        splog, 'Writing '+jpegfile+'.jpg'
        nw_rgb_make, iim_nw, rim_nw, gim_nw, name=jpegfile, $
          scales=scales, nonlinearity=3.0, rebinfactor=2, quality=95
    endif
    return
end
