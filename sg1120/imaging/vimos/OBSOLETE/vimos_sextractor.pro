pro vimos_sextractor, dec03=dec03, feb06=feb06
; jm07jan24nyu - run SExtractor on all the images in preparation for
;                SCAMP and SWARP

    if ((not keyword_set(dec03)) and (not keyword_set(feb06))) or $
      (keyword_set(dec03) and keyword_set(feb06)) then begin
       splog, 'Must specify either DEC03 or FEB06 keyword!'
       return
    endif
    
    sexpath = sg1120_path(/sex)
    datapath = vimos_path(dec03=dec03,feb06=feb06)+'sg1120/'

    imagelist = file_search(datapath+'ra.sg1120*_[B,V,R].fits',count=nimage)
;   imagelist = file_search(datapath+'ra.sg1120*_*_Q?.fits',count=nimage)
    weightlist = repstr(imagelist,'.fits','.weight.fits')
    catlist = repstr(imagelist,'.fits','.cat')

    sexconfig = sexpath+'sg1120.sex'
    sexparam = sexpath+'sg1120.sex.param'
    sexconv = sexpath+'default.conv'
    sexnnw = sexpath+'default.nnw'

    t0 = systime(1)
    for ii = 0L, nimage-1L do begin

       if (file_test(imagelist[ii],/regular) eq 0L) then begin
          splog, 'Image '+imagelist+ 'not found.'
          return
       endif

       print, 'SExtracting '+imagelist[ii]

       hdr = headfits(imagelist[ii],ext=1)
       seeing = string(sxpar(hdr,'SEEING'),format='(G0.0)') ; approximate value for the observing run
       if (seeing le 0.0) then seeing = '1.0' ; temporary fix!
       gain = string(sxpar(hdr,'GAIN'),format='(G0.0)')

       spawn, 'sex '+imagelist[ii]+' -c '+sexconfig+' -CATALOG_NAME '+catlist[ii]+' -CATALOG_TYPE FITS_LDAC'+$
         ' -PARAMETERS_NAME '+sexparam+' -FILTER_NAME '+sexconv+' -STARNNW_NAME '+sexnnw+$
         ' -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE '+weightlist[ii]+' -WEIGHT_THRESH 0'+$
         ' -GAIN '+gain+' -SEEING_FWHM '+seeing+' -INTERP_TYPE NONE'+$
         ' -MEMORY_BUFSIZE 8192 -VERBOSE_TYPE NORMAL', /sh

    endfor
    splog, 'Total time to run = ', (systime(1)-t0)/60.0, ' minutes.'
    
return
end
    
