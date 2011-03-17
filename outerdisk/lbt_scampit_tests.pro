pro lbt_scampit, sextractor=sextractor, scamp=scamp, swarp=swarp, test=test, $
  striptnx=striptnx, jpeg=jpeg, astrometry=astrometry, update_astrometry=update_astrometry
; jm07jan25nyu - generate all the mosaics

; scamp p94_u_10.ext4.cat -CDSCLIENT_EXEC aclient -CHECKPLOT_DEV NULL -MOSAIC_TYPE SHARE_PROJAXIS -SAVE_REFCATALOG Y -MERGEDOUTCAT_TYPE FITS_LDAC -ASTREF_CATALOG SDSS-R3
    
    rootpath = '/global/bias1/ioannis/lbt/rawdata/'
;   datapath = '/global/bias1/ioannis/lbt/data_4ext/'
    datapath = '/global/bias1/ioannis/lbt/data/'
    sexpath = '/global/bias1/ioannis/lbt/sex/'
    mosaicpath = '/global/bias1/ioannis/lbt/'
    astrom_path = datapath+'astrometry.net/'
    
    scampconfig = sexpath+'mosaic.scamp'
    swarpconfig = sexpath+'mosaic.swarp'
    sexconfig = sexpath+'mosaic.sex'
    sexparam = sexpath+'mosaic.sex.param'
    sexconv = sexpath+'default.conv'
    sexnnw = sexpath+'default.nnw'

; ---------------------------------------------------------------------------    

    if keyword_set(striptnx) then begin
       splog, 'Recommend you run UPDATE_ASTROMETRY after this!'
       imagelist = [file_search(rootpath+'p94_*_?.fits'),file_search(rootpath+'p94_*_10.fits')]
       newimagelist = datapath+file_basename(imagelist)
       t0 = systime(1)
;      for ii = 4L, 4L do begin
       for ii = 0L, n_elements(imagelist)-1L do begin
          bighdr = headfits(imagelist[ii],ext=0)
;         mwrfits, 0, newimagelist[ii], /create ; NOTE!
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
             sxdelpar, hdr, 'COMMENT'
             sxdelpar, hdr, 'DATE'
             sxaddpar, hdr, 'DATE_OBS', sxpar(bighdr,'DATE_OBS')
             sxaddpar, hdr, 'OBJECT', sxpar(oldhdr,'OBJECT')
             sxaddpar, hdr, 'GAIN', float(sxpar(oldhdr,'GAIN'))
             sxaddpar, hdr, 'RDNOISE', float(sxpar(oldhdr,'RDNOISE'))
             sxaddpar, hdr, 'EXPTIME', float(sxpar(oldhdr,'EXPTIME'))
             sxaddpar, hdr, 'AIRMASS', float(sxpar(bighdr,'AIRMASS'))
             sxaddpar, hdr, 'FILTER', sxpar(bighdr,'FILTER')
             sxdelpar, hdr, 'HISTORY'
             putast, hdr, a
; sky-subtract
             splog, '   Sky-subtracting...'
             case iext of
                1L: imsky = im[0:500,*]
                2L: imsky = im[*,0:500]
                3L: imsky = im[nx-501:nx-1,*]
                4L: imsky = im[*,ny-501:ny-1]
             endcase
             mmm, imsky, skymode, skysig;, /debug
             imnosky = im - skymode
; reject cosmic rays
             splog, '   Rejecting cosmic rays...'
             imsig = dsigma(im)
             invvar = 0.0*im+1.0/(skysig^2.0)
             reject_cr, imnosky, invvar, [0.496,0.246], rejects, $
               nrejects=nrejects, c2fudge=c2fudge, niter=10L
             splog, '   Identified '+string(nrejects,format='(I0)')+' cosmic rays.'
             invvar[rejects] = 0.0
             if (nrejects gt 0L) then imnosky = djs_maskinterp(imnosky,(invvar le 0),iaxis=0,/const)
; write out             
             outfile = repstr(newimagelist[ii],'.fits','.ext'+string(iext,format='(I0)')+'.fits')
             splog, 'Writing '+outfile
             mwrfits, float(imnosky), outfile, hdr, /create
;            splog, 'Writing '+newimagelist[ii]+', extension '+string(iext,format='(I0)')
;            mwrfits, im, newimagelist[ii], hdr
;            if (iext eq 4L) then begin
;               splog, 'Writing '+newimagelist_ext4[ii]
;               mwrfits, im, newimagelist_ext4[ii], hdr, /create
;            endif else begin
;               splog, 'Writing '+newimagelist[ii]+', extension '+string(iext,format='(I0)')
;               mwrfits, im, newimagelist[ii], hdr
;            endelse
          endfor
       endfor 
       splog, 'Total time to pre-process images = ', (systime(1)-t0)/60.0, ' minutes.'
    endif

; ---------------------------------------------------------------------------    

    if keyword_set(astrometry) then begin
       splog, 'Remove all FITS files from '+astrom_path+' [Y/N]?'
       cc = get_kbrd(1)
       if (strupcase(cc) eq 'Y') then begin
          flist = [file_search(astrom_path+'*.fits'),file_search(astrom_path+'*.html')]
          rmfile, flist
       endif
       imagelist = file_search(datapath+'p94_*_*.ext?.fits')
       pushd, astrom_path
       spawn, 'hansel.csh -e john.moustakas@gmail.com --units arcmin --field-width 7.67 '+$
         '--pixel-scale 0.00373 --uncertainty 10 '+strjoin(imagelist,' '), /sh
       popd
    endif
    
; ---------------------------------------------------------------------------    

    if keyword_set(update_astrometry) then begin
;      imagelist = file_search(datapath+'p94_*_*.ext?.fits')
       imagelist = file_search(datapath+'p94_u_*.ext?.fits')
       wcslist = file_search(astrom_path+repstr(file_basename(imagelist),'.fits','*_wcs.fits'))
       if (n_elements(imagelist) ne n_elements(wcslist)) then begin
          splog, 'IMAGELIST and WCSLIST must match.'
          return
       endif else begin
          niceprint, imagelist, wcslist
          splog, 'Does this look OK?'
          cc = get_kbrd(1)
          if (strupcase(cc) ne 'Y') then return
       endelse
       t0 = systime(1)
       for ii = 0L, n_elements(imagelist)-1L do begin
          splog, 'Updating header for '+imagelist[ii]
          newhdr = headfits(wcslist[ii])
          sxdelpar, newhdr, 'COMMENT'
          modfits, imagelist[ii], 0, newhdr
       endfor
       splog, 'Total time to update headers = ', (systime(1)-t0)/60.0, ' minutes.'
    endif

; ---------------------------------------------------------------------------    

    if keyword_set(sextractor) then begin

       imagelist = file_search(datapath+'p94_u_*.ext?.fits')
       weightlist = repstr(imagelist,'.fits','.weight.fits')
       catlist = repstr(imagelist,'.fits','.cat')
       seglist = repstr(imagelist,'.fits','.seg.fits')

       t0 = systime(1)
;      for ii = 0L, 2L do begin
       for ii = 0L, n_elements(imagelist)-1L do begin
          print, 'SExtracting '+imagelist[ii]
          spawn, 'sex '+imagelist[ii]+' -c '+sexconfig+' -CATALOG_NAME '+catlist[ii]+' -CATALOG_TYPE FITS_LDAC'+$
            ' -PARAMETERS_NAME '+sexparam+' -DETECT_THRESH 1.5 -ANALYSIS_THRESH 1.5 -FILTER_NAME '+sexconv+' -STARNNW_NAME '+sexnnw+$
            ' -WEIGHT_TYPE NONE -WEIGHT_IMAGE '+weightlist[ii]+' -WEIGHT_THRESH 0'+$ ; NOTE! weight map
            ' -VERBOSE_TYPE NORMAL -NTHREADS 4'+$
            ' -CHECKIMAGE_TYPE NONE -CHECKIMAGE_NAME '+seglist[ii], /sh
;           ' -CHECKIMAGE_TYPE SEGMENTATION -CHECKIMAGE_NAME '+seglist[ii], /sh
       endfor
       splog, 'Total time to generate SE catalogs = ', (systime(1)-t0)/60.0, ' minutes.'
    endif

; ---------------------------------------------------------------------------    

    if keyword_set(scamp) then begin
    
       catlist = file_search(datapath+'p94_u_10.ext?.cat')
;      catlist = file_search(datapath+'p94_u_10.ext4.cat')

       astref_catalog = 'USNO-B1'
       astrinstru_key = 'FILTER'
       photinstru_key = 'FILTER'
       magzero_key = 'PHOT_C'
       extinct_key = 'PHOT_K'
       photomflag_key = 'PHOTFLAG'

       t0 = systime(1)
       for iter = 0L, 2L do begin

          case iter of
             0L: begin
                degree = '1'
                mosaic_type = 'UNCHANGED'
;               mosaic_type = 'SHARE_PROJAXIS'
                position_maxerr = '1.0'
                posangle_maxerr = '3.0'
                aheader_suffix = '.ahead'
             end
             1L: begin
                degree = '2'
                mosaic_type = 'UNCHANGED' ; 'SHARE_PROJAXIS'
                position_maxerr = '0.2'
                posangle_maxerr = '3.0'
                aheader_suffix = '.head'
             end
             else: begin
                degree = '3'
                mosaic_type = 'UNCHANGED' ; 'FIX_FOCALPLANE' ; 'SHARE_PROJAXIS'
                position_maxerr = '0.1'
                posangle_maxerr = '0.5'
                aheader_suffix = '.head'
             end
          endcase

          spawn, 'scamp '+strjoin(catlist,',')+' -c '+scampconfig+' -CHECKPLOT_DEV PSC -CDSCLIENT_EXEC aclient'+$
            ' -PIXSCALE_MAXERR 1.1 -POSANGLE_MAXERR '+posangle_maxerr+' -POSITION_MAXERR '+position_maxerr+$
            ' -ASTREF_CATALOG SDSS-R3'+$
;           ' -AHEADER_GLOBAL '+datapath+'lbt.ahead '+$
            ' -AHEADER_SUFFIX '+aheader_suffix+' -DISTORT_DEGREES '+degree+' -MOSAIC_TYPE '+mosaic_type+$
            ' -ASTRINSTRU_KEY '+astrinstru_key+' -NTHREADS 3 -WRITE_XML N -VERBOSE_TYPE NORMAL', /sh

;         spawn, 'scamp '+strjoin(catlist,',')+' -c '+scampconfig+' -CHECKPLOT_DEV NULL -CDSCLIENT_EXEC aclient'+$
;           ' -ASTREF_CATALOG '+astref_catalog+' -SAVE_REFCATALOG N -MERGEDOUTCAT_TYPE NONE'+$ 
;           ' -PIXSCALE_MAXERR 1.1 -POSANGLE_MAXERR '+posangle_maxerr+' -POSITION_MAXERR '+position_maxerr+' -MOSAIC_TYPE '+mosaic_type+$
;           ' -ASTRINSTRU_KEY '+astrinstru_key+' -DISTORT_DEGREES '+degree+$
;           ' -SOLVE_PHOTOM Y -PHOTINSTRU_KEY '+photinstru_key+' -MAGZERO_KEY '+magzero_key+$
;           ' -EXTINCT_KEY '+extinct_key+' -PHOTOMFLAG_KEY '+photomflag_key+$
;           ' -WRITE_XML N -VERBOSE_TYPE NORMAL -NTHREADS 2 -AHEADER_SUFFIX '+aheader_suffix, /sh ; NOTE!
stop
       endfor 
       splog, 'Total time to run scamp = ', (systime(1)-t0)/60.0, ' minutes.'

    endif
       
; ---------------------------------------------------------------------------    
    
    if keyword_set(swarp) then begin

; global parameters       

       if keyword_set(test) then begin
          radec = '00:00:00.0,+00:00:00.0'
          center_type = 'ALL' 
          pixelscale_type = 'MANUAL' 
          pixel_scale = '0.224'
          image_size = '0'
          header_only = 'Y'
          suffix = 'test'
       endif else begin
          radec = '12:50:53.0,+41:07:14.0' ; '12:51:00.0,+41:04:40.0'
          center_type = 'MANUAL' 
          pixelscale_type = 'MANUAL' 
          pixel_scale = '0.224'
          image_size = '8500,7590' ; '7000,5000' ; '7850,6150'
          header_only = 'N'
          suffix = ''
       endelse

       combinetype = 'MEDIAN' ; 'WEIGHTED' ; 'MEDIAN' ; 'WEIGHTED'
       gain_keyword = 'GAIN'
       projection_err = '0.001'
       keywords = 'OBJECT,FILTER'

       bandpass = ['U','V']
       mosaic_file = mosaicpath+'M94_'+bandpass+suffix+'.fits'
       mosaic_weightfile = repstr(mosaic_file,'.fits','.weight.fits')
       jpegfile = mosaicpath+'M94_UV.jpeg'

       for ib = 0L, 0L do begin
;      for ib = 1L, 1L do begin
;      for ib = 0L, 1L do begin

          if file_test(mosaic_file[ib],/regular) then begin
             splog, 'Mosaic '+mosaic_file[ib]+' exists.'
;            return
          endif
          
          case ib of
             0L: imagelist = file_search(datapath+'p94_u_7.fits')
;            0L: imagelist = [file_search(datapath+'p94_u_?.fits'),file_search(datapath+'p94_u_10.fits')]
             1L: imagelist = [file_search(datapath+'p94_v_?.fits'),file_search(datapath+'p94_v_10.fits')]
             else:
          endcase

          t0 = systime(1)
          splog, 'Building '+mosaic_file[ib]+' from '+string(n_elements(imagelist),format='(I0)')+' images.'
          spawn, 'swarp '+strjoin(imagelist,',')+' -c '+swarpconfig+' -HEADER_ONLY '+header_only+$
            ' -IMAGEOUT_NAME '+mosaic_file[ib]+' -WEIGHTOUT_NAME '+mosaic_weightfile[ib]+$
            ' -WEIGHT_TYPE NONE -WEIGHT_SUFFIX _weight.fits -WEIGHT_THRESH 0 -HEADER_SUFFIX .head'+$
            ' -COMBINE_TYPE '+combinetype+' -BLANK_BADPIXELS Y -CELESTIAL_TYPE EQUATORIAL'+$
            ' -PROJECTION_TYPE TAN -PROJECTION_ERR '+projection_err+' -CENTER_TYPE '+center_type+' -CENTER '+radec+$
            ' -PIXELSCALE_TYPE '+pixelscale_type+' -PIXEL_SCALE '+pixel_scale+' -IMAGE_SIZE '+image_size+$
            ' -RESAMPLING_TYPE LANCZOS3 -INTERPOLATE Y -GAIN_KEYWORD '+gain_keyword+' -SUBTRACT_BACK N'+$
;           ' -RESAMPLE_DIR resamp -DELETE_TMPFILES N'+$
            ' -COMBINE_BUFSIZE 1024 -COPY_KEYWORDS '+keywords+' -WRITE_FILEINFO N -VERBOSE_TYPE NORMAL -NTHREADS 2', /sh
          splog, 'Total time to build '+mosaic_file[ib]+' = ', (systime(1)-t0)/60.0, ' minutes.'

stop          
          
       endfor
          
    endif

; ---------------------------------------------------------------------------    

    if keyword_set(jpeg) then begin

;      splog, 'Reading '+mosaicpath+'M94_U.fits'
;      image = mrdfits(mosaicpath+'M94_U.fits',0,hdr,/silent)
;
;      pagemaker, nx=1, ny=1, xspace=0.0, yspace=0.0, xmargin=1.0, ymargin=1.0, $
;        /normal, position=pos
;
;      imsize = size(image,/dimension)
;      xsize = imsize[0] & xcen = xsize/2.0 & ysize = imsize[1] & ycen = ysize/2.0
;      xaxis = (findgen(xsize)-xcen) ;*pixscale ; [arcsec]
;      yaxis = (findgen(ysize)-ycen) ;*pixscale ; [arcsec]
;
;      img = logscl(image,negative=1,exp=0.5)
;      img = asinhscl(image,negative=1,alpha=5.0,beta=16.0,omin=0,omax=250)
;      
;      plotimage, img, /normal, position=pos, margin=0, imgxrange=minmax(xaxis), $
;        imgyrange=minmax(yaxis), xtickname=replicate(' ',10), $
;        ytickname=replicate(' ',10), xsty=5, ysty=5;, /preserve
;
;          
;      splog, 'Reading '+mosaicpath+'M94_V.fits'
;      Vim = mrdfits(mosaicpath+'M94_V.fits',0,hdr,/silent)
       
       jpegfile = mosaicpath+'M94_UV.jpeg'
       splog, 'Reading '+mosaicpath+'M94_U.fits'
       Uim = mrdfits(mosaicpath+'M94_U.fits',0,hdr,/silent)
       imsize = size(Uim,/dim)
       nx = 512 * (min(imsize)/2048)
       Uim_nw = (temporary(Uim))$
         [imsize[0]/2L-nx/2L:imsize[0]/2L+nx/2L-1L, $
         imsize[1]/2L-nx/2L:imsize[1]/2L+nx/2L-1L]
       splog, 'Reading '+mosaicpath+'M94_V.fits'
       Vim_nw = (mrdfits(mosaicpath+'M94_V.fits',0,/silent))$
         [imsize[0]/2L-nx/2L:imsize[0]/2L+nx/2L-1L, $
         imsize[1]/2L-nx/2L:imsize[1]/2L+nx/2L-1L]
       splog, 'Writing '+jpegfile

       nw_rgb_make, Vim_nw, Vim_nw, Uim_nw, name=jpegfile, $
         scales=[2.0,1.0,1.0], nonlinearity=3.0, rebinfactor=2, quality=95
       
stop       
       
    endif
   
return       
end
