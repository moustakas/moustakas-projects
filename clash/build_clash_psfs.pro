pro clash_psf_sex, imfile, catalog_name=catalog_name, zpt=zpt, gain=gain
; build a basic SE catalog using the F160W band in double-image mode 
;   rmsfile = repstr(imfile,'drz','rms')
    weightfile = repstr(imfile,'drz','wht')

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
    
pro build_clash_psfs, sexcatalogs=sexcatalogs
; jm13may30siena - build the PSFs for CLASH in all 16 bands using all
; the clusters 

    archive = getenv('CLASH_ARCHIVE')
    clash = rsex(getenv('CLASH_DIR')+'/clash_sample.sex')
    clash = clash[8] ; use just MACS1206 for now
    ncl = n_elements(clash)

    filt = bcgimf_filterlist(short=short,weff=weff,instr=instr)
    reffilt = (where(short eq 'f160w'))[0] ; reference filter
    nfilt = n_elements(filt)

    mosaicpath = getenv('IM_ARCHIVE_DIR')+'/'+cluster+'/images/'
;   mosaicpath = getenv('IM_ARCHIVE_DIR')+'/'+cluster+'/patched_images/'
    catpath = getenv('IM_ARCHIVE_DIR')+'/'+cluster+'/catalogs/'
    propath = getenv('CLASH_DIR')+'/streams/'
    streamspath = getenv('IM_PROJECTS_DIR')+'/clash/streams/'
    psfpath = streamspath+'psf/'

    parameters_file = propath+'sex.param'
    



stop
    
    magpivot = 25.0
    int = 27.0
    slope = 0.9
    magfaint = 25.0
    magaxis = range(15,magfaint,500)
    for ic = 0, ncl-1 do begin
       if file_test(archive+strtrim(clash[ic].dirname,2),/dir) then begin
          splog, 'Working on cluster '+clash[ic].cluster

          catpath = archive+strtrim(clash[ic].dirname,2)+'/HST/catalogs/'+$
            'mosaicdrizzle_image_pipeline/IR_detection/SExtractor/'
          mosaicpath = archive+strtrim(clash[ic].dirname,2)+'/HST/images/'+$
            'mosaicdrizzle_image_pipeline/scale_65mas/'
          cat = rsex(catpath+'detectionImage.cat')
          cat.mag_auto += 25
          cat.mag_aper += 25

          istar = where(cat.mag_auto gt 0.0 and cat.mag_auto lt magfaint and $
            cat.mag_aper lt poly(cat.mag_auto-magpivot,[int,slope]) and $
            cat.flags ne 0,nstar)
          plot, [0], [0], /nodata, xrange=[15,30], yrange=[35,15], $
            xtitle='mag_auto (AB mag)', ytitle='mag_aper (AB mag)'
          djs_oplot, cat.mag_auto, cat.mag_aper, psym=8
          djs_oplot, cat[istar].mag_auto, cat[istar].mag_aper, psym=8, color='green'
          djs_oplot, magaxis, poly(magaxis-magpivot,[int,slope]), line=0
          djs_oplot, magfaint*[1,1], [poly(magfaint-magpivot,[int,slope]),!y.crange[1]]
;         legend, strupcase(short[ii]), /right, /top, box=0

stop          
          
; get cutouts
          imfile = clash[ic].shortname+'_mosaic_065mas_wfc3ir_f160w_drz_????????.fits.gz'
          image = mrdfits(mosaicpath+imfile,0,hdr)
          extast, hdr, astr
          ad2xy, cat[istar].alpha_j2000, cat[istar].delta_j2000, astr, xx, yy
          for is = 0, nstar-1 do begin
             hextract, image, hdr, cutout, newhdr, xx[is]-20, $
               xx[is]+20, yy[is]-20, yy[is]+20, /silent
             atv, cutout, /bl
          endfor
          
          
stop          

;         for ii = 0, nfilt-1 do begin
;            cat = rsex(catpath+strtrim(clash[ic].shortname,2)+'_'+short[ii]+'.cat')
;            istar = where(cat.mag_auto gt 0.0 and cat.mag_auto lt magfaint and $
;              cat.mag_aper lt poly(cat.mag_auto-magpivot,[int,slope]))
;            plot, [0], [0], /nodata, xrange=[15,30], yrange=[35,15], $
;              xtitle='mag_auto (AB mag)', ytitle='mag_aper (AB mag)'
;            djs_oplot, cat.mag_auto, cat.mag_aper, psym=8
;            djs_oplot, cat[istar].mag_auto, cat[istar].mag_aper, psym=8, color='green'
;
;            djs_oplot, magaxis, poly(magaxis-magpivot,[int,slope]), line=0
;            djs_oplot, magfaint*[1,1], [poly(magfaint-magpivot,[int,slope]),!y.crange[1]]
;
;            legend, strupcase(short[ii]), /right, /top, box=0
;         endfor
       endif else splog, 'No data for cluster '+clash[ic].cluster
    endfor

return
end
    
