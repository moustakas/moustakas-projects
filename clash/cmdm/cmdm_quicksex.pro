pro cmdm_quicksex, imagefile, weightfile=weightfile, $
  write_regfile=write_regfile, silent=silent
; jm12feb28ucsd - build a quick SE catalog

    sexpath = clash_path()+'projects/cmdm/data/'

    paramfile = sexpath+'quicksex.param'
    convfile = sexpath+'default.conv'
    nnwfile = sexpath+'default.nnw'

    cl = 'macs1206'
    imfile = sexpath+'macsj1206_2009_IC_O1.fits'
    weightfile = sexpath+'macsj1206_2009_IC_O1.weight.fits'
    catfile = repstr(imfile,'.fits','.cat')
    
    nimage = n_elements(imfile)
    config = init_sex_config(nimage)

    config.catalog_name = catfile
    config.parameters_name = paramfile
    config.filter_name = convfile
    config.starnnw_name = nnwfile

    config.phot_fluxfrac = ['0.5,0.9,0.95']

    config.catalog_type = 'ASCII_HEAD'
;   config.catalog_type = 'FITS_LDAC'
    config.detect_thresh = 1.5
    config.analysis_thresh = 1.5

    config.weight_type = 'MAP_WEIGHT'
    config.weight_image = weightfile
    
    config.seeing_fwhm = 0.7
    config.mag_zeropoint = sxpar(headfits(imfile),'MAGZPT')
    config.checkimage_type = 'NONE'    ; SEGMENTATION
;   config.checkimage_name = seglist

    im_sex, imfile, config, silent=silent ; do not pass CONFIGFILE

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
    
return
end
