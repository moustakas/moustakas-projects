pro mips_sextractor, hizzle=hizzle
; jm07jan24nyu - run SExtractor on all the images in preparation for
;   SCAMP and SWARP

    sexpath = sg1120_path(/sex)
    datapath = mips_path()

    if keyword_set(hizzle) then begin
    
; first crop the MIPS mosaic to be the same as the Ks-band mosaic     

       im = mrdfits('mips_mosaic_01292007.fits',0,h,/silent)
       mwrfits, im*0.0+1.0, 'mips_mosaic_01292007.weight.fits', h, /create
       
       extast, headfits(flamingos_path(/mosaics)+'sg1120_ks_weighted_order5.fits'), mosaic_astr
       hizzle, 'mips_mosaic_01292007.fits', outfile='mips_mosaic_01292007_ksmatched.fits', $
         weight_suffix='.weight', mosaic_astr=mosaic_astr, degree=2

       return
       
    endif

    mipsimage = datapath+'mips_mosaic_01292007_ksmatched.fits'
    mipsweightimage = datapath+'mips_mosaic_01292007_ksmatched.weight.fits'
    ksimage = flamingos_path(/mosaics)+'sg1120_ks_weighted_order5.fits'
    
    catlist = repstr(mipsimage,'.fits','.cat')

    sexconfig = sexpath+'sg1120.sex'
    sexparam = sexpath+'sg1120.sex.param'
    sexconv = sexpath+'default.conv'
    sexnnw = sexpath+'default.nnw'

    print, 'SExtracting '+mipsimage
    spawn, 'sex '+ksimage+' '+mipsimage+' -c '+sexconfig+' -CATALOG_NAME '+catlist+' -CATALOG_TYPE FITS_LDAC'+$
      ' -PARAMETERS_NAME '+sexparam+' -DETECT_THRESH 1.5 -FILTER_NAME '+sexconv+' -STARNNW_NAME '+sexnnw+$
      ' -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE '+mipsweightimage+' -WEIGHT_THRESH 0'+$ ; NOTE! weight map
      ' -GAIN 1.0 -PIXEL_SCALE 0.32 -SEEING_FWHM 6.0 -BACK_TYPE MANUAL -BACKPHOTO_TYPE GLOBAL -INTERP_TYPE NONE'+$
      ' -MEMORY_BUFSIZE 8192 -VERBOSE_TYPE NORMAL -NTHREADS 2', /sh

return
end
    
