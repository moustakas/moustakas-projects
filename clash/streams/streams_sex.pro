pro dostreams_sex, imfile, catalog_name=catalog_name, zpt=zpt, $
  gain=gain, psfname=psfname, parameters_file=parameters_file
; build the SExtractor catalog; many of the parameters are adopted
; from D. Coe's choice of parameters
    rmsfile = repstr(imfile,'drz','rms')
    weightfile = repstr(imfile,'drz','wht')

    config = init_sex_config()
    config.catalog_name = catalog_name
    config.psf_name = psfname

    config.parameters_name = parameters_file
    config.filter_name = '/usr/local/share/sextractor-2.8.6/config/gauss_2.0_5x5.conv'
    config.starnnw_name = getenv('IMPRO_DIR')+'/etc/default.nnw'

    config.clean_param = 1.2

    config.analysis_thresh = 2.0
    config.detect_thresh = 2.0
    config.deblend_mincont = 0.0015
    config.deblend_nthresh = 32
    config.detect_minarea = 9

    config.checkimage_type = 'SEGMENTATION,BACKGROUND'
    config.checkimage_name = strjoin(repstr(catalog_name,'.cat')+'_'+['segm','back']+'.fits',',')

    config.phot_apertures = 15

; background subtraction    
    config.backphoto_type = 'GLOBAL'
    config.backphoto_thick = 24
    config.back_filtersize = 5
    config.back_size = 128
    
    config.catalog_type = 'ASCII_HEAD' ; 'FITS_LDAC'
;   config.weight_type = 'MAP_RMS'
;   config.weight_image = rmsfile
    config.weight_type = 'MAP_WEIGHT'
    config.weight_image = weightfile
    
    config.mag_zeropoint = zpt
    config.satur_level = 249500
    config.gain = gain
    config.seeing_fwhm = 0.140
;   config.verbose_type = 'FULL'

    im_sex, imfile, config, silent=silent ; do not pass CONFIGFILE

return
end
    
pro dostreams_sex_psfex, imfile, catalog_name=catalog_name, zpt=zpt, gain=gain
; build the pre-PSFEx SExtractor catalog
    rmsfile = repstr(imfile,'drz','rms')
;   weightfile = repstr(imfile,'drz','wht')

    config = init_sex_config()
    config.catalog_name = catalog_name
    config.parameters_name = getenv('IMPRO_DIR')+'/etc/psfex.sex.param'
    config.filter_name = getenv('IMPRO_DIR')+'/etc/psfex.sex.conv'
    config.starnnw_name = getenv('IMPRO_DIR')+'/etc/default.nnw'
    
    config.detect_minarea = 4
    config.detect_thresh = 4
    config.analysis_thresh = 4
    config.phot_apertures = 15
    
    config.catalog_type = 'FITS_LDAC'
    config.weight_type = 'MAP_RMS'
    config.weight_image = rmsfile
;   config.weight_type = 'MAP_WEIGHT'
;   config.weight_image = weightfile
    
    config.mag_zeropoint = zpt
    config.satur_level = 249500
    config.gain = gain

    im_sex, imfile, config, silent=silent ; do not pass CONFIGFILE

return
end
    
pro dostreams_psfex, sexcatalog, psfpath=psfpath
; run PSFEx
    config = init_psfex_config()

;   config.checkimage_type = 'RESIDUALS'
;   config.checkimage_name = psfpath+'resi.fits'
    config.checkimage_type = 'CHI,PROTOTYPES,SAMPLES,RESIDUALS,SNAPSHOTS'
    config.checkimage_name = strjoin(psfpath+['chi.fits','proto.fits','samp.fits','resi.fits','snap.fits'],',')

    config.basis_type = 'NONE'
    config.sample_variability = 0.1
    config.badpixel_filter = 'Y'
    
    config.psfvar_nsnap = 5
    config.psfvar_degrees = 2
;   config.psfvar_degrees = 0 ; constant PSF
    config.psf_dir = psfpath
;   config.psfvar_keys = ''   ; no spatial dependence
;   config.sample_fwhmrange = '0.1,3' ; expected range of stellar FWHM
;   config.verbose_type = 'FULL'
    im_psfex, sexcatalog, config, silent=silent
return
end
    
pro streams_sex, get_psf=get_psf, build_catalogs=build_catalogs
; jm13may23siena - use SExtractor to fit MACS1206

    cluster = 'macs1206'
    mosaicpath = getenv('IM_ARCHIVE_DIR')+'/'+cluster+'/images/'
;   mosaicpath = getenv('IM_ARCHIVE_DIR')+'/'+cluster+'/patched_images/'
    catpath = getenv('IM_ARCHIVE_DIR')+'/'+cluster+'/catalogs/'
    propath = getenv('CLASH_DIR')+'/streams/'
    streamspath = getenv('IM_PROJECTS_DIR')+'/clash/streams/'
    psfpath = streamspath+'psf/'

    parameters_file = propath+'sex.param'
    
; cluster and filter properties    
    clash = rsex(getenv('CLASH_DIR')+'/clash_sample.sex')
    this =  where('macs1206' eq strtrim(clash.cluster_short,2))

    filt = bcgimf_filterlist(short=short,instr=instr,weff=weff,zpt=zpt)
    nfilt = n_elements(filt)

; sort through each filter
    for ib = 12, nfilt-1 do begin
;   for ib = 0, nfilt-1 do begin

; some global filenames
       catalog_name = streamspath+cluster+'_'+short[ib]+'.cat'
       psfname = psfpath+cluster+'_'+short[ib]+'.psf'
       
; for GAIN see notes at
; http://www.ifa.hawaii.edu/~rgal/science/sextractor_notes.html 
       imfile = mosaicpath+'macs1206_mosaic_065mas_'+instr[ib]+$
         '_'+short[ib]+'_drz_20110815.fits'
       hdr = headfits(imfile)
;      gain = sxpar(hdr,'exptime')        
       gain = sxpar(hdr,'ccdgain')*sxpar(hdr,'exptime')        

; make a test image       
       im = mrdfits(imfile,0,hdr)
       rms = mrdfits(repstr(imfile,'drz','rms'),0,rmshdr)
       hextract, im, hdr, testim, testhdr, 2300, 3200, 2900, 3400
       hextract, rms, rmshdr, testrms, testrmshdr, 2300, 3200, 2900, 3400
       mwrfits, testim, streamspath+'test_drz_f160w.fits', testhdr, /create
       mwrfits, testrms, streamspath+'test_rms_f160w.fits', testrmshdr, /create
       testimfile = streamspath+'test_drz_f160w.fits'
       
; --------------------------------------------------
; build the pre-PSFEx SExtractor catalog
       if keyword_set(get_psf) then begin
          psfcatalog_name = psfpath+cluster+'_'+short[ib]+'.cat'
          dostreams_sex_psfex, imfile, catalog_name=psfcatalog_name, $
            zpt=zpt[ib], gain=gain
          dostreams_psfex, psfcatalog_name, psfpath=psfpath
       endif

; --------------------------------------------------
; build the SE catalog       
       if keyword_set(build_catalogs) then begin
;         dostreams_sex, testimfile, catalog_name=catalog_name, $
          dostreams_sex, imfile, catalog_name=catalog_name, $
            psfname=psfname, zpt=zpt[ib], gain=gain, $
            parameters_file=parameters_file
          cc = rsex(catalog_name)
       endif
stop       
       
; optionally write a region file
       if keyword_set(write_regfile) then begin
          for ii = 0, n_elements(catfile)-1 do begin
             for jj = 1, 36 do begin
                cat1 = mrdfits(catfile[ii],2*jj,/silent)
                if (jj eq 1) then cat = cat1 else cat = [cat,cat1]
             endfor
             write_ds9_regionfile, cat.xwin_world, cat.ywin_world, color='red', $
               file=repstr(catfile[ii],'.cat','.reg')
          endfor
       endif

    endfor

return
end
    
