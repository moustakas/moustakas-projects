pro merge_redmapper_isedfit, centrals=centrals, satellites=satellites
; jm14sep16siena - merge the iSEDfit and redmapper catalogs, split
;   into satellites and BCGs  

    catpath = redmapper_path(ver=ver,/catalogs)
    isedpath = redmapper_path(/isedfit)
    datapath = redmapper_path()

; pick out the centrals & satellites; note that PHOT contains
; everything in the members catalog, as well as my photometry 
    phot = mrdfits(catpath+'redmapper_'+ver+'_phot.fits.gz',1)
    icen = where(phot.isbcg eq 1)
    isat = where(phot.isbcg eq 0)

    ised = mrdfits(isedpath+'redmapper_fsps_v2.4_miles_'+$
      'chab_charlot_sfhgrid01.fits.gz',1)
    kcorr = mrdfits(catpath+'redmapper_'+ver+'_kcorr.fits.gz',1)

; for the centrals need to sort the BCGs file
    bcgs = mrdfits(catpath+'dr8_run_redmapper_'+ver+'_lgt5_catalog.fits.gz',1)

    match, phot[icen].mem_match_id, bcgs.mem_match_id, m1, m2
    srt = sort(m1)
    bcgs = bcgs[m2[srt]]
;   print, minmax(bcgs.ra-mems[icen].ra), minmax(bcgs.dec-mems[icen].dec)
;   help, where(bcgs.z_lambda-mems[icen].z ne 0)

    centrals = struct_addtags(struct_trimtags(phot[icen],$
      select=['mem_match_id','z','ra','dec','r',$
      'p','pfree','theta_i','theta_r']),$
      struct_trimtags(bcgs,except=['ra','dec',$
      'mem_match_id','model_mag','model_magerr','photoid','z']))
    centrals = struct_addtags(temporary(centrals),$
      struct_trimtags(kcorr[icen],except=['z','maggies','ivarmaggies','bestmaggies']))
    centrals = struct_addtags(temporary(centrals),$
      struct_trimtags(ised[icen],except=['z','ra','dec']))

    im_mwrfits, centrals, datapath+'redmapper_isedfit_'+ver+'_centrals.fits', /clobber

; now do the satellites    
    satellites = struct_trimtags(phot[isat],except=['maggies','ivarmaggies',$
      'model_mag','model_magerr','photoid','isbcg'])
    satellites = struct_addtags(temporary(satellites),struct_trimtags(kcorr[isat],$
      except=['z','maggies','ivarmaggies','bestmaggies']))
    satellites = struct_addtags(temporary(satellites),$
      struct_trimtags(ised[isat],except=['z','ra','dec']))

    im_mwrfits, satellites, datapath+'redmapper_isedfit_'+ver+'_satellites.fits', /clobber

return
end
