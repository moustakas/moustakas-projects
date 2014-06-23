pro objectmask_sex, imfile, sexpath=sexpath, catalog_name=catalog_name, $
  weightfile=weightfile, detect_imfile=detect_imfile, zpt=zpt, gain=gain
; build a basic SE catalog using the F160W band in double-image mode 

    path = file_dirname(imfile)+'/'
    checkfile = repstr(imfile,'.fits','-segm.fits')

    config = init_sex_config()
    config.catalog_name = catalog_name
    config.parameters_name = sexpath+'objectmask_sex.param'
    config.filter_name = sexpath+'default.conv'
    config.starnnw_name = sexpath+'default.nnw'
    
    config.detect_minarea = 5.0
    config.detect_thresh = 10
    config.analysis_thresh = 10
;   config.phot_apertures = 1
    
    config.checkimage_type = 'SEGMENTATION'
    config.checkimage_name = checkfile

    config.catalog_type = 'FITS_LDAC'
;   config.weight_type = 'MAP_RMS'
;   config.weight_image = rmsfile
;   config.weight_type = 'NONE'
    config.weight_type = 'MAP_WEIGHT'
    config.weight_image = weightfile
    
;   config.mag_zeropoint = zpt
    config.satur_level = 249500
    config.pixel_scale = 0.065
    config.gain = gain

    im_sex, imfile, config;, detect_imagelist=detect_imfile

return
end
    
pro bcgmstar_objectmask, copyfiles=copyfiles, sextractor=sextractor, makemask=makemask
; jm14may21siena - build an object mask for the BCGs/Mstar project  

; note! images are converted to [10^-12 erg/s/cm^2/Hz] (pico-maggies)
    sample = read_bcgmstar_sample()
    ncl = n_elements(sample) 

; psf path (see BUILD_CLASH_PSFS)
    rootpath = bcgmstar_path(/objectmask)
    psfpath = getenv('IM_ARCHIVE_DIR')+'/projects/clash/psfs/'
    
; specifiy the filters and some other handy info    
    filt = bcgmstar_filterlist(short=short,instr=instr,$
      weff=weff,zpt=zpt)
    allfiltinfo = replicate({filt: '', short: '', instr: '', $
      weff: 0.0, zpt: 0.0},n_elements(filt))
    allfiltinfo.filt = filt
    allfiltinfo.short = short
    allfiltinfo.instr = instr
    allfiltinfo.weff = weff
    allfiltinfo.zpt = zpt

    splog, 'HACK!!!'
    allfiltinfo = allfiltinfo[13]
    nfilt = n_elements(allfiltinfo)
    struct_print, allfiltinfo

; wrap on each cluster    
;   for ic = ncl-1, ncl-1 do begin
    for ic = 0, ncl-1 do begin
       cluster = strtrim(sample[ic].shortname,2)
       splog, 'Working on cluster '+cluster
       outpath = rootpath+cluster+'/'
       if file_test(outpath,/dir) eq 0 then file_mkdir, outpath

;; read the PSF stars in the reference band; choosing a masking radius
;; that varies with magnitude
;; faint - 12 pixels
;; 17th magnitude - 60 pixels
;; 19th magnitude - 35 pixels
;       stars = mrdfits(psfpath+cluster+'/'+cluster+'-f160w-stars.fits.gz',1)
;       starrad = poly(stars.mag_auto,[308,-14.0])>15.0
;       nstar = n_elements(stars)

;; mask the stars
;          adxy, hdr, stars.ra, stars.dec, xx, yy
;          mx1 = weighted_quantile(imagecube[*,*,ib],quant=0.9)
;          cgimage, imagecube[*,*,ib], clip=3, /negative, stretch=5, minvalue=0.0, maxvalue=mx1, $
;            margin=0, /keep_aspect, /save
;          for ii = 0, nstar-1 do tvcircle, starrad[ii], xx[ii], yy[ii], color='cyan', /data
       
; --------------------------------------------------
; copy the files we need
       if keyword_set(copyfiles) then begin
          splog, 'Copying files...'

; not all clusters have all filters, so check for that case here
          bcgmodelpath = getenv('CLASH_ARCHIVE')+'/'+strtrim(sample[ic].dirname,2)+$
            '/HST/galaxy_subtracted_images/marc/'
          mosaicpath = getenv('CLASH_ARCHIVE')+'/'+strtrim(sample[ic].dirname,2)+$
            '/HST/images/mosaicdrizzle_image_pipeline/scale_65mas/'
          these = where(file_test(mosaicpath+cluster+'_mosaic_065mas_'+$
            strtrim(allfiltinfo.instr,2)+'_'+strtrim(allfiltinfo.short,2)+'_drz_????????.fits*'),nfilt)
          if nfilt gt n_elements(allfiltinfo) or nfilt eq 0 then message, 'Problem here!'
          filtinfo = allfiltinfo[these]
          reffilt = where(filtinfo.short eq 'f160w') ; reference filter

; read the original mosaics and inverse variance maps in a datacube 
          drzfiles = file_search(mosaicpath+cluster+'_mosaic_065mas_'+$
            filtinfo.instr+'_'+filtinfo.short+'_drz_????????.fits*')
          whtfiles = file_search(mosaicpath+cluster+'_mosaic_065mas_'+$
            filtinfo.instr+'_'+filtinfo.short+'_wht_????????.fits*')

;         bcgmodelfiles = file_search(bcgmodelpath+cluster+$
;           '_mosaic_065mas_*_'+filtinfo.short+'_drz_*_model.fits.gz')
          bcgmodelfiles = file_search(bcgmodelpath+cluster+$
            '_mosaic_065mas_*_'+filtinfo.short+'_drz_*_BCG.fits.gz')
;         bcgmodelfiles = file_search(bcgmodelpath+cluster+$
;           '_mosaic_065mas_*_'+filtinfo.short+'_drz_*_subt.fits.gz')

          for ib = nfilt-1, 0, -1 do begin
             spawn, 'gunzip -c '+drzfiles[ib]+' > '+$
               outpath+cluster+'-'+filtinfo[ib].short+'-image.fits', /sh
             spawn, 'gunzip -c '+whtfiles[ib]+' > '+$
               outpath+cluster+'-'+filtinfo[ib].short+'-ivar.fits', /sh
;            spawn, 'gunzip -c '+bcgmodelfiles[ib]+' > '+$
;              outpath+cluster+'-nobcg-'+filtinfo[ib].short+'.fits', /sh
             mwrfits, mrdfits(drzfiles[ib],0,hdr,/silent)-$
               mrdfits(bcgmodelfiles[ib],0,/silent), outpath+cluster+'-'+filtinfo[ib].short+$
               '-nobcg.fits', hdr, /create
          endfor
       endif

; --------------------------------------------------
; create the sextractor segmentation maps
       if keyword_set(sextractor) then begin
          splog, 'Building segmentation images...'
          imfile = file_search(outpath+cluster+'-*-nobcg.fits',count=nfilt)
          weightfile = repstr(imfile,'-nobcg','-ivar')
          catalog_name = repstr(imfile,'.fits','.cat')
          for ib = 0, nfilt-1 do begin
             hdr = headfits(imfile[ib])
             gain = sxpar(hdr,'ccdgain')*sxpar(hdr,'exptime')        
             objectmask_sex, imfile[ib], catalog_name=catalog_name[ib], $
               weightfile=weightfile, sexpath=rootpath, gain=gain
          endfor
       endif 
       
; --------------------------------------------------
; build the final object mask
       if keyword_set(makemask) then begin
          splog, 'Building the object mask...'
          catfile = file_search(outpath+cluster+'-*-nobcg.cat',count=nfilt)
          segmfile = repstr(catfile,'.cat','-segm.fits')
          imfile = repstr(catfile,'-nobcg.cat','-image.fits')
          maskfile = repstr(catfile,'-nobcg.cat','-objectmask.fits')
          testfile = repstr(catfile,'-nobcg.cat','-masked.fits')

          for ib = 0, nfilt-1 do begin
             segm = mrdfits(segmfile,0,hdr,/silent)
             bigmask = segm*0
             sz = size(segm,/dim)

; sometimes the residuals in the BCG subtraction get detected as
; sources, which messes up our ellipse-fitting; fix that here
             case cluster of
                'a209': segm[sz[0]/2-40:sz[0]/2+40,sz[1]/2-40:sz[1]/2+40] = 0
                'a383': segm[where(segm eq 290)] = 0 ; residuals
                'a1423': segm[where(segm eq 155)] = 0 ; galaxy very close to the core!
                'macs1206': segm[where(segm eq 282 or segm eq 285)] = 0 ; residuals
; these are probably all real sources (esp 216), but the masking is
; too aggressive
                'clj1226': segm[where(segm eq 232 or segm eq 231 or segm eq 216)] = 0 
                'macs1720': segm[where(segm eq 290 or segm eq 305)] = 0 
                'a2261': segm[where(segm eq 160 or segm eq 196 or $
                  segm eq 222 or segm eq 219)] = 0 ; residuals
; these are probably all real sources, but the masking is too
; aggressive 
;               'macs2129': segm[where(segm eq 250 or segm eq 252)] = 0
                'macs2129': segm[where(segm eq 252)] = 0 
                'rxj2129': segm[where(segm eq 218 or segm eq 222 or segm eq 220)] = 0 
                'ms2137': segm[sz[0]/2-20:sz[0]/2+20,sz[1]/2-20:sz[1]/2+20] = 0
; 201 is a real galaxy, but it's smack in the center!
                'rxj2248': segm[where(segm eq 201)] = 0 
                else:
             endcase
             
             cat = mrdfits(catfile[ib],2)
             big = where(cat.isoarea_image gt 1000,nbig,comp=small,ncomp=nsmall)
             splog, nbig, nsmall
             if nbig eq 0 or nsmall eq 0 then message, 'Problem!'
;            small = where(cat.isoarea_image lt 200,nsmall)
;            if nsmall eq 0 then message, 'Problem!'
;            for ii = 0L, nsmall-1 do segm[where(segm eq cat[small[ii]].number)] = 0
             for ii = 0L, nbig-1 do bigmask[where(segm eq cat[big[ii]].number)] = 1

             bigmask = dsmooth(bigmask,10)
;            bigmask = dsmooth(dsmooth(bigmask,10) gt 0,10) gt 0
;            bigmask = dsmooth(bigmask,10) gt 0
             mask = (dsmooth(segm,2) + bigmask) gt 0
;            mask = (dsmooth(segm,5) + bigmask) gt 0
             splog, 'Writing '+maskfile
             mwrfits, mask, maskfile, hdr, /create

; write out the masked original image so we can check it
             testim = mrdfits(imfile,0,/silent)*(mask eq 0)
             mwrfits, testim, testfile, hdr, /create 
stop             
          endfor 
       endif
    endfor                      ; close cluster loop

return
end


;             hdr = (mrdfits(catfile[ib],1)).field_header_card
;             cat = mrdfits(catfile[ib],2)
;             cat = cat[reverse(sort(cat.a_image))]
;             sz = [sxpar(hdr,'NAXIS1'),sxpar(hdr,'NAXIS2')]
;             mask = intarr(sz)
;             ximage = findgen(sz[0])#(1+fltarr(sz[1]))
;             yimage = (1+fltarr(sz[0]))#findgen(sz[1])
;             for ii = 0, 20 do begin
;;            for ii = 0, n_elements(cat)-1 do begin
;                print, ii
;                indx1 = get_ellipse_indices(reform(ximage,cmproduct(sz)),$
;                  reform(yimage,cmproduct(sz)),major=cat[ii].a_image,$
;                  minor=cat[ii].b_image,angle=cat[ii].theta_image+90,$
;                  xcenter=cat[ii].x_image, ycenter=cat[ii].y_image)
;;               dist_ellipse, mask1, sz, cat[ii].x_image, cat[ii].y_image, $
;;                 cat[ii].elongation, cat[ii].theta_image+90
;;               mask += (mask1 lt 3*cat[ii].a_image) ge 1
;                mask[indx1] = 1
;             endfor
;             mask = mask ge 1
;             mwrfits, mask, 'junk.fits', hdr, /create
