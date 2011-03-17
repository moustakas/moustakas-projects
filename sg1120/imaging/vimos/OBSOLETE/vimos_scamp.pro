pro vimos_scamp, usno=usno, dec03=dec03, feb06=feb06
; jm07jan24nyu - run SCAMP

    if ((not keyword_set(dec03)) and (not keyword_set(feb06))) or $
      (keyword_set(dec03) and keyword_set(feb06)) then begin
       splog, 'Must specify either DEC03 or FEB06 keyword!'
       return
    endif
    
    sexpath = sg1120_path(/sex)
    datapath = vimos_path(dec03=dec03,feb06=feb06)+'sg1120/'

;   catlist = file_search(datapath+'ra.sg1120*_B.cat')
    catlist = file_search(datapath+'ra.sg1120*_[B,V,R].cat')
;   catlist = file_search(datapath+'ra.sg1120*_*_Q?.cat')

    root = 'vimos'
    scampconfig = sexpath+'sg1120.scamp'

    if keyword_set(usno) then astref_catalog = 'USNO-B1' else astref_catalog = '2MASS'
    refout_catpath = datapath
    mergedoutcat_name = datapath+root+'.scamp.cat.fits'
    astrinstru_key = 'FILTER'
    photinstru_key = 'FILTER'
    magzero_key = 'PHOT_C'
    extinct_key = 'PHOT_K'
    photomflag_key = 'PHOTFLAG'
    degree = '3'

    for iter = 0L, 2L do begin

       if (iter eq 0L) then begin
          mosaic_type = 'UNCHANGED'
          position_maxerr = '3.0'
          aheader_suffix = '.ahead'
       endif else begin
          mosaic_type = 'FIX_FOCALPLANE'
          position_maxerr = '0.5'
          aheader_suffix = '.head'
       endelse

       spawn, 'scamp '+strjoin(catlist,',')+' -c '+scampconfig+$
         ' -ASTREF_CATALOG '+astref_catalog+' -SAVE_REFCATALOG Y -REFOUT_CATPATH '+refout_catpath+$
         ' -MERGEDOUTCAT_NAME '+mergedoutcat_name+' -MERGEDOUTCAT_TYPE FITS_LDAC '+$
         ' -PIXSCALE_MAXERR 1.2 -POSANGLE_MAXERR 2.0 -POSITION_MAXERR '+position_maxerr+' -MOSAIC_TYPE '+mosaic_type+$
         ' -ASTRINSTRU_KEY '+astrinstru_key+' -DISTORT_DEGREES '+degree+$
         ' -SOLVE_PHOTOM Y -PHOTINSTRU_KEY '+photinstru_key+' -MAGZERO_KEY '+magzero_key+$
         ' -EXTINCT_KEY '+extinct_key+' -PHOTOMFLAG_KEY '+photomflag_key+$
         ' -WRITE_XML N -VERBOSE_TYPE NORMAL -AHEADER_SUFFIX '+aheader_suffix, /sh ; NOTE!

    endfor    
    
return
end
