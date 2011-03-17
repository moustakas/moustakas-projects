pro vimos_swarp, dec03=dec03, feb06=feb06, test=test
; jm07jan24nyu - run SWARP; note: the (ra,dec) are hard-coded here; to
;   allow SWARP to choose an arbitrary center you have to change
;   CENTER_TYPE to ALL

    if ((not keyword_set(dec03)) and (not keyword_set(feb06))) or $
      (keyword_set(dec03) and keyword_set(feb06)) then begin
       splog, 'Must specify either DEC03 or FEB06 keyword!'
       return
    endif
    
    sexpath = sg1120_path(/sex)
    mosaicpath = vimos_path(/mosaics)
    datapath = vimos_path(dec03=dec03,feb06=feb06)+'sg1120/'
    if keyword_set(dec03) then suffix = '2003' else suffix = '2006'
    
    swarpconfig = sexpath+'sg1120.swarp'

    if keyword_set(test) then begin
       radec = '00:00:00.0,+00:00:00.0'
       center_type = 'ALL' 
       pixelscale_type = 'MANUAL' 
       pixel_scale = strtrim(0.205,2)
       image_size = '0'
       header_only = 'Y'
       jpeg = 0L
    endif else begin
       if strmatch(suffix,'*2003*') then begin
          radec = '11:20:03.0,-12:04:30.0'
          center_type = 'MANUAL' 
          pixelscale_type = 'MANUAL' 
          pixel_scale = strtrim(0.205,2)
          image_size = '5842,5120'
       endif else begin
          radec = '11:19:54.0,-12:02:43.0'
          center_type = 'MANUAL' 
          pixelscale_type = 'MANUAL' 
          pixel_scale = strtrim(0.205,2)
          image_size = '5833,5133'
       endelse
       header_only = 'N'
       jpeg = 1L
    endelse

    combinetype = 'WEIGHTED'    ; 'MEDIAN' ; 'WEIGHTED'
    gain_keyword = 'GAIN'
    projection_err = '0.001'
    keywords = 'OBJECT,FILTER'
    
; B-, V-, and R-band mosaics    

    bandpass = ['B','V','R']
    mosaic_file = mosaicpath+'sg1120_'+bandpass+'_'+suffix+'.fits'
;   mosaic_file = mosaicpath+'sg1120_'+bandpass+'_'+strlowcase(combinetype)+'_'+suffix+'.fits'
    mosaic_weightfile = repstr(mosaic_file,'.fits','.weight.fits')
    jpegfile = mosaicpath+'sg1120_bvr_'+suffix+'.jpeg'
;   jpegfile = mosaicpath+'sg1120_bvr_'+strlowcase(combinetype)+'_'+suffix+'.jpeg'

    for ib = 0L, n_elements(bandpass)-1L do begin

; run SWARP; generate the image file list based on what was used in
; SCAMP
       
       catheadlist = file_search(datapath+'ra.sg1120*_'+bandpass[ib]+'.head')
;      catheadlist = file_search(datapath+'ra.sg1120*_'+bandpass[ib]+'_Q?.head')

       imagelist = repstr(catheadlist,'.head','.fits')
       weightlist = repstr(imagelist,'.fits','.weight.fits')

       splog, 'Building '+mosaic_file[ib]+' from '+string(n_elements(imagelist),format='(I0)')+' images.'
       spawn, 'swarp '+strjoin(imagelist,',')+' -c '+swarpconfig+' -HEADER_ONLY '+header_only+$
         ' -IMAGEOUT_NAME '+mosaic_file[ib]+' -WEIGHTOUT_NAME '+mosaic_weightfile[ib]+$
         ' -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_SUFFIX .weight.fits -WEIGHT_THRESH 0'+$
         ' -COMBINE_TYPE '+combinetype+' -BLANK_BADPIXELS Y -CELESTIAL_TYPE EQUATORIAL'+$
         ' -PROJECTION_TYPE TAN -PROJECTION_ERR '+projection_err+' -CENTER_TYPE '+center_type+' -CENTER '+radec+$
         ' -PIXELSCALE_TYPE '+pixelscale_type+' -PIXEL_SCALE '+pixel_scale+' -IMAGE_SIZE '+image_size+$
         ' -RESAMPLING_TYPE LANCZOS3 -INTERPOLATE N -GAIN_KEYWORD '+gain_keyword+' -SUBTRACT_BACK Y'+$
         ' -COMBINE_BUFSIZE 512 -COPY_KEYWORDS '+keywords+' -WRITE_FILEINFO N -VERBOSE_TYPE NORMAL -NTHREADS 2', /sh

    endfor

; generate the color mosaic

    if keyword_set(jpeg) then begin

       bim = mrdfits(mosaic_file[0],0,bhdr,/silent)
       vim = mrdfits(mosaic_file[1],0,vhdr,/silent)
       rim = mrdfits(mosaic_file[2],0,rhdr,/silent)
       gim = (bim + rim) / 2.0
       imsize = size(bim,/dim)

       nx = 512 * (min(imsize)/512)

       bim_nw = bim[imsize[0]/2L-nx/2L:imsize[0]/2L+nx/2L-1L,imsize[1]/2L-nx/2L:imsize[1]/2L+nx/2L-1L]
       vim_nw = vim[imsize[0]/2L-nx/2L:imsize[0]/2L+nx/2L-1L,imsize[1]/2L-nx/2L:imsize[1]/2L+nx/2L-1L]
       rim_nw = rim[imsize[0]/2L-nx/2L:imsize[0]/2L+nx/2L-1L,imsize[1]/2L-nx/2L:imsize[1]/2L+nx/2L-1L]
       gim_nw = gim[imsize[0]/2L-nx/2L:imsize[0]/2L+nx/2L-1L,imsize[1]/2L-nx/2L:imsize[1]/2L+nx/2L-1L]

       splog, 'Writing '+jpegfile
       djs_rgb_make, rim_nw, gim_nw, bim_nw, name=jpegfile, scales=[1.0,1.0,1.0], $
         nonlinearity=3.0, rebinfactor=2, quality=75

    endif

return
end
    
