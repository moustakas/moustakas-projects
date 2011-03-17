pro lbt_mosaics, preproc=preproc, sextractor=sextractor, scamp=scamp, swarp=swarp, test=test, jpeg=jpeg
; J. Moustakas, 2007 May - generate mosaics from Dennis' LBT imaging of M94

    rootpath = '/global/bias1/ioannis/lbt/'
    rawpath = rootpath+'/rawdata/'
    datapath = rootpath+'data/'
    sexpath = rootpath+'sex/'
    mosaicpath = rootpath
    
    scampconfig = sexpath+'default.scamp'
    swarpconfig = sexpath+'default.swarp'
    sexconfig = sexpath+'default.sex'
    sexparam = sexpath+'lbt.sex.param'
    sexconv = sexpath+'default.conv'
    sexnnw = sexpath+'default.nnw'

; ---------------------------------------------------------------------------    
; pre-process the mosaics (much of this should be incorporated in the
; preliminary reductions): (1) re-generate the astrometry headers,
; which are *not FITS compliant*; (2) generate a simple bad pixel map
; (basically mask out the bad columns, which otherwise messes with SE;
; (3) do a simple sky-subtraction; (4) reject cosmic rays; (5) write out
; ---------------------------------------------------------------------------    
    
    if keyword_set(preproc) then begin
       imagelist = [file_search(rawpath+'p94_*_?.fits'),file_search(rawpath+'p94_*_10.fits')]
       newimagelist = datapath+file_basename(imagelist)
       weightlist = repstr(newimagelist,'.fits','.weight.fits')
       t0 = systime(1)
       for ii = 0L, n_elements(imagelist)-1L do begin
          bighdr = headfits(imagelist[ii],ext=0)
          mwrfits, 0, newimagelist[ii], /create
          mwrfits, 0, weightlist[ii], /create
          for iext = 1L, 4L do begin
             splog, 'Processing '+imagelist[ii]+', extension '+string(iext,format='(I0)')
             im = mrdfits(imagelist[ii],iext,oldhdr,/silent)
             imsize = size(im,/dim) & nx = imsize[0] & ny = imsize[1]
             orientation = -(sxpar(oldhdr,'CROTA1') MOD 360.0)
; HACK: all images get same pixscale
             if (not keyword_set(pixscale)) then pixscale = abs(sxpar(oldhdr,'CDELT1'))*3600.0D ; [arcsec/pixel]
             a = hogg_make_astr(sxpar(oldhdr,'CRVAL1'),sxpar(oldhdr,'CRVAL2'),$
               sxpar(oldhdr,'NAXIS1')*pixscale/3600.0D,sxpar(oldhdr,'NAXIS2')*pixscale/3600.0D,$
               orientation=orientation,pixscale=pixscale/3600.0D)
             a.crpix = [sxpar(oldhdr,'CRPIX1'),sxpar(oldhdr,'CRPIX2')] ; NOTE!
             mkhdr, hdr, im ; generate a basic FITS header
             gain = float(sxpar(oldhdr,'GAIN'))
             rdnoise = float(sxpar(oldhdr,'RDNOISE'))
             sxdelpar, hdr, 'COMMENT'
             sxdelpar, hdr, 'DATE'
             sxaddpar, hdr, 'DATE_OBS', sxpar(bighdr,'DATE_OBS')
             sxaddpar, hdr, 'OBJECT', sxpar(oldhdr,'OBJECT')
             sxaddpar, hdr, 'GAIN', gain
             sxaddpar, hdr, 'RDNOISE', rdnoise
             sxaddpar, hdr, 'EXPTIME', float(sxpar(oldhdr,'EXPTIME'))
             sxaddpar, hdr, 'AIRMASS', float(sxpar(bighdr,'AIRMASS'))
             sxaddpar, hdr, 'FILTER', sxpar(bighdr,'FILTER')
             sxdelpar, hdr, 'HISTORY'
             putast, hdr, a
; generate a weight map
             splog, '   Generating a weight map...'
             badpix = im*0.0 ; 0=good, 1=bad
             case iext of
                1L: begin
                   badpix[61,*] = 1.0
                   badpix[375,*] = 1.0
                   badpix[581,*] = 1.0
                   badpix[702,*] = 1.0
                   badpix[1050,*] = 1.0
                   badpix[1067,*] = 1.0
                   badpix[1384,*] = 1.0
                   badpix[1648,*] = 1.0
                   badpix[*,4605:4607] = 1.0
                end
                2L: begin
                   badpix[0,99:4607] = 1.0
                   badpix[1827:1830,99:4607] = 1.0
                   badpix[2047,99:4607] = 1.0
                   badpix[*,4605:4607] = 1.0
                end
                3L: begin
                   badpix[0,*] = 1.0
                   badpix[2047,99:4607] = 1.0
                   badpix[1618:1621,3712:4607] = 1.0
                   badpix[1405,4361] = 1.0
                   badpix[*,4605:4607] = 1.0
                end
                4L: begin
                   badpix[0,*] = 1.0
                   badpix[2047,99:4607] = 1.0
                   badpix[21:29,2496] = 1.0
                   badpix[42:43,2471] = 1.0
                   badpix[*,4605:4607] = 1.0
                end
             endcase
             badpix = smooth(badpix,5,/nan,/edge) ; grow the mask
             invmap = badpix eq 0.0 ; 1=good, 0=bad
;; sky-subtract
             splog, '   Sky-subtracting...'
             imsky = im[1700:1900,100:200] ; from Dennis
             mmm, imsky, skymode, skysig ;, /debug
             imnosky = im - skymode
;            imnosky = im
; reject cosmic rays
             splog, '   Rejecting cosmic rays...'
             imsig = dsigma(imnosky)
             invvar = 0.0*imnosky+1.0/(imsig^2.0)
             reject_cr, imnosky, invvar, [0.496,0.246], rejects, $
               nrejects=nrejects, c2fudge=c2fudge, niter=10L
             splog, '   Identified '+string(nrejects,format='(I0)')+' cosmic rays.'
             if (nrejects gt 0L) then begin
                invmap[rejects] = 0.0
                imnosky = imnosky*(invmap gt 0.0)
;               imnosky = djs_maskinterp(imnosky,(invmap le 0),iaxis=0,/const)
             endif
; write out             
             splog, 'Writing '+newimagelist[ii]+', extension '+string(iext,format='(I0)')
             mwrfits, float(imnosky), newimagelist[ii], hdr
             splog, 'Writing '+weightlist[ii]+', extension '+string(iext,format='(I0)')
             mwrfits, byte(invmap), weightlist[ii], hdr
          endfor
       endfor 
       splog, 'Total time to pre-process images = ', (systime(1)-t0)/60.0, ' minutes.'
    endif

; ---------------------------------------------------------------------------    
; generate SE catalogs
; ---------------------------------------------------------------------------    
    
    if keyword_set(sextractor) then begin

       imagelist = [file_search(datapath+'p94_*_?.fits'),file_search(datapath+'p94_*_10.fits')]
       weightlist = repstr(imagelist,'.fits','.weight.fits')
       catlist = repstr(imagelist,'.fits','.cat')
       seglist = repstr(imagelist,'.fits','.seg.fits')

       t0 = systime(1)
       for ii = 0L, n_elements(imagelist)-1L do begin
          print, 'SExtracting '+imagelist[ii]
          spawn, 'sex '+imagelist[ii]+' -c '+sexconfig+' -CATALOG_NAME '+catlist[ii]+' -CATALOG_TYPE FITS_LDAC'+$
            ' -PARAMETERS_NAME '+sexparam+' -DETECT_THRESH 2.0 -ANALYSIS_THRESH 2.0 -FILTER_NAME '+sexconv+' -STARNNW_NAME '+sexnnw+$
            ' -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE '+weightlist[ii]+' -WEIGHT_THRESH 0'+$
            ' -VERBOSE_TYPE NORMAL -NTHREADS 4 -DEBLEND_MINCONT 0.01 -BACK_TYPE AUTO -BACK_VALUE 0.0'+$
            ' -MEMORY_BUFSIZE 4096 -CHECKIMAGE_TYPE NONE -CHECKIMAGE_NAME '+seglist[ii], /sh
;           ' -CHECKIMAGE_TYPE SEGMENTATION -CHECKIMAGE_NAME '+seglist[ii], /sh
       endfor
       splog, 'Total time to generate SE catalogs = ', (systime(1)-t0)/60.0, ' minutes.'
    endif

; ---------------------------------------------------------------------------    
; do the astrometry, using SCAMP
; ---------------------------------------------------------------------------    
    
    if keyword_set(scamp) then begin
    
       catlist = file_search(datapath+'p94_*_*.cat')

       astref_catalog = 'USNO-B1' & astref_band = 'Rf'
;      astref_catalog = 'SDSS-R5' & astref_band = 'r'
       astrinstru_key = 'FILTER'
       photinstru_key = 'FILTER'
       magzero_key = 'PHOT_C'
       extinct_key = 'PHOT_K'
       photomflag_key = 'PHOTFLAG'

       t0 = systime(1)
       for iter = 0L, 2L do begin ; three iterations is good

          case iter of
             0L: begin
                degree = '1'
                mosaic_type = 'LOOSE'
;               mosaic_type = 'UNCHANGED'
                position_maxerr = '3.0'
                posangle_maxerr = '2.0'
                aheader_suffix = '.ahead'
             end
             1L: begin
                degree = '3'
;               mosaic_type = 'LOOSE'
                mosaic_type = 'FIX_FOCALPLANE'
                position_maxerr = '0.5'
                posangle_maxerr = '1.0'
                aheader_suffix = '.head'
             end
             else: begin
                degree = '3'
                mosaic_type = 'FIX_FOCALPLANE'
                position_maxerr = '0.5'
                posangle_maxerr = '1.0'
                aheader_suffix = '.head'
             end
          endcase

          spawn, 'scamp '+strjoin(catlist,' ')+' -c '+scampconfig+' -CHECKPLOT_DEV PSC'+$
            ' -PIXSCALE_MAXERR 1.1 -POSANGLE_MAXERR '+posangle_maxerr+' -POSITION_MAXERR '+position_maxerr+$
            ' -ASTREF_CATALOG '+astref_catalog+' -ASTREF_BAND '+astref_band+' -MERGEDOUTCAT_TYPE FITS_LDAC -SAVE_REFCATALOG Y'+$
            ' -AHEADER_SUFFIX '+aheader_suffix+' -DISTORT_DEGREES '+degree+' -MOSAIC_TYPE '+mosaic_type+$
            ' -CROSSID_RADIUS 5.0 -FWHM_THRESHOLDS 3.0,100.0 -STABILITY_TYPE EXPOSURE'+$
            ' -ASTRINSTRU_KEY '+astrinstru_key+' -NTHREADS 4 -WRITE_XML N -VERBOSE_TYPE NORMAL', /sh

       endfor 
       splog, 'Total time to run scamp = ', (systime(1)-t0)/60.0, ' minutes.'

    endif
       
; ---------------------------------------------------------------------------    
; finally generate the big mosaics; use identical astrometric headers
; for both the U- and V-band
; ---------------------------------------------------------------------------    
    
    if keyword_set(swarp) then begin

       if keyword_set(test) then begin
          radec = '00:00:00.0,+00:00:00.0'
          center_type = 'ALL' 
          pixelscale_type = 'MANUAL' 
          pixel_scale = '0.224'
          image_size = '0'
          header_only = 'Y'
          suffix = 'test'
       endif else begin
          radec = '12:50:52.0,+41:08:30.0'
;         radec = '12:50:53.0,+41:07:14.0' ; M94/NED coordinates
          center_type = 'MANUAL' 
          pixelscale_type = 'MANUAL' 
          pixel_scale = '0.224'
          image_size = '7900,8300'
          header_only = 'N'
          suffix = ''
       endelse

;      combinetype = 'AVERAGE'
       combinetype = 'MEDIAN' ; better to use AVERAGE
       gain_keyword = 'GAIN'
       projection_err = '0.001'
       keywords = 'OBJECT,FILTER' ; keywords to copy

       bandpass = ['U','V']
       mosaic_file = mosaicpath+'M94_'+bandpass+suffix+'.fits'
       mosaic_weightfile = repstr(mosaic_file,'.fits','.weight.fits')

       for ib = 0L, 1L do begin

          if file_test(mosaic_file[ib],/regular) then begin
             splog, 'WARNING! Mosaic '+mosaic_file[ib]+' exists -- press Y to overwrite:'
             cc = get_kbrd(1)
             if (strupcase(cc) ne 'Y') then return 
          endif
          
          case ib of
             0L: imagelist = [file_search(datapath+'p94_u_?.fits'),file_search(datapath+'p94_u_10.fits')]
             1L: imagelist = [file_search(datapath+'p94_v_?.fits'),file_search(datapath+'p94_v_10.fits')]
             else:
          endcase

          t0 = systime(1)
          splog, 'Building '+mosaic_file[ib]+' from '+string(n_elements(imagelist),format='(I0)')+' images.'
          spawn, 'swarp '+strjoin(imagelist,',')+' -c '+swarpconfig+' -HEADER_ONLY '+header_only+$
            ' -IMAGEOUT_NAME '+mosaic_file[ib]+' -WEIGHTOUT_NAME '+mosaic_weightfile[ib]+$
            ' -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_SUFFIX .weight.fits -WEIGHT_THRESH 0 -HEADER_SUFFIX .head'+$
            ' -COMBINE_TYPE '+combinetype+' -BLANK_BADPIXELS Y -CELESTIAL_TYPE EQUATORIAL'+$
            ' -PROJECTION_TYPE TAN -PROJECTION_ERR '+projection_err+' -CENTER_TYPE '+center_type+' -CENTER '+radec+$
            ' -PIXELSCALE_TYPE '+pixelscale_type+' -PIXEL_SCALE '+pixel_scale+' -IMAGE_SIZE '+image_size+$
            ' -RESAMPLING_TYPE LANCZOS3 -INTERPOLATE Y -GAIN_KEYWORD '+gain_keyword+' -SUBTRACT_BACK N'+$
;           ' -RESAMPLE_DIR resamp -DELETE_TMPFILES N'+$
            ' -COMBINE_BUFSIZE 1024 -COPY_KEYWORDS '+keywords+' -WRITE_FILEINFO N -VERBOSE_TYPE NORMAL -NTHREADS 4', /sh
          splog, 'Total time to build '+mosaic_file[ib]+' = ', (systime(1)-t0)/60.0, ' minutes.'

       endfor
          
    endif

; ---------------------------------------------------------------------------    
; test code to generate a color mosaic; unfortunately, we don't have a
; red bandpass, and the V-band image is saturated, so it doesn't work
; very well
; ---------------------------------------------------------------------------    

    if keyword_set(jpeg) then begin

       jpegfile = mosaicpath+'M94_UV.jpeg'
       splog, 'Reading '+mosaicpath+'M94_U.fits'
       Uim = mrdfits(mosaicpath+'M94_U.fits',0,hdr,/silent)
       imsize = size(Uim,/dim)
;      xcen = imsize/2L & ycen = imsize/2L
       xcen = 3930L & ycen = 3830L
       nx = 512 * (min(imsize)/2048)
       Uim_nw = (temporary(Uim))[xcen-nx/2L:xcen+nx/2L-1L,ycen-nx/2L:ycen+nx/2L-1L]
       splog, 'Reading '+mosaicpath+'M94_V.fits'
       Vim_nw = (mrdfits(mosaicpath+'M94_V.fits',0,/silent))$
         [xcen-nx/2L:xcen+nx/2L-1L,ycen-nx/2L:ycen+nx/2L-1L]
       splog, 'Writing '+jpegfile

       nw_rgb_make, Vim_nw, Vim_nw, Uim_nw, name=jpegfile, $
         scales=[1.0,1.0,5.0], nonlinearity=3.0, rebinfactor=2, quality=95
       
stop       
       
    endif
   
return       
end
