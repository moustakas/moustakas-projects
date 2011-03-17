pro mips_scamp, usno=usno
; jm07jan24nyu - run SCAMP

    sexpath = sg1120_path(/sex)
    datapath = mips_path()
;   catlist = file_search(datapath+'mips_mosaic_01292007.cat',count=nimage)
    catlist = [datapath+'mips_mosaic_01292007_ksmatched.cat',$
      flamingos_path(/mosaics)+'sg1120_ks_weighted_order5.cat']
    
    root = 'mips'
    scampconfig = sexpath+'sg1120.scamp'

    if keyword_set(usno) then astref_catalog = 'USNO-B1' else astref_catalog = '2MASS'
    refout_catpath = datapath
    mergedoutcat_name = root+'.scamp.cat'
    mosaic_type = 'UNCHANGED'
    astrinstru_key = '1,2'
    photinstru_key = '1,2'
;   astrinstru_key = 'FILTER'
;   photinstru_key = 'FILTER'
    magzero_key = 'PHOT_C'
    extinct_key = 'PHOT_K'
    photomflag_key = 'PHOTFLAG'
    degree = '2'

    for iter = 0L, 1L do begin 

       if (iter eq 0L) then aheader_suffix = '.ahead' else aheader_suffix = '.head'
       case iter of
          0L: position_maxerr = '0.5'
          1L: position_maxerr = '0.5'
          else: position_maxerr = '0.5'
       endcase

       spawn, 'scamp '+strjoin(catlist,',')+' -c '+scampconfig+$
         ' -ASTREF_CATALOG '+astref_catalog+' -SAVE_REFCATALOG Y -REFOUT_CATPATH '+refout_catpath+$
         ' -MERGEDOUTCAT_NAME '+mergedoutcat_name+' -MERGEDOUTCAT_TYPE FITS_LDAC'+$
         ' -PIXSCALE_MAXERR 1.2 -POSANGLE_MAXERR 2.0 -POSITION_MAXERR '+position_maxerr+' -MOSAIC_TYPE '+mosaic_type+$
         ' -CROSSID_RADIUS 4 -MATCH_FLIPPED N '+$
         ' -ASTRINSTRU_KEY '+astrinstru_key+' -DISTORT_DEGREES '+degree+' -ASTRCLIP_NSIGMA 5.0'+$
         ' -SOLVE_PHOTOM N -PHOTINSTRU_KEY '+photinstru_key+' -MAGZERO_KEY '+magzero_key+$
         ' -EXTINCT_KEY '+extinct_key+' -PHOTOMFLAG_KEY '+photomflag_key+$
         ' -WRITE_XML N -VERBOSE_TYPE NORMAL -AHEADER_SUFFIX '+aheader_suffix, /sh ; NOTE!

stop
       
    endfor
       
return
end
