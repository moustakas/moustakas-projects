pro read_redmapper, centrals=centrals, satellites=satellites
; jm13aug27siena - read all the redmapper outputs

    datapath = redmapper_path(ver=ver)

; pick out the centrals & satellites; note that PHOT contains
; everything in the members catalog, as well as my photometry 
    if n_elements(centrals) eq 0L or n_elements(satellites) eq 0L then begin
       phot = mrdfits(datapath+'redmapper_'+ver+'_phot.fits.gz',1)
       icen = where(phot.isbcg eq 1)
       isat = where(phot.isbcg eq 0)

       ised = mrdfits(datapath+'redmapper_fsps_v2.4_miles_'+$
         'chab_charlot_sfhgrid01.fits.gz',1)

       satellites = struct_addtags(struct_trimtags(phot[isat],$
         except=['maggies','ivarmaggies']),struct_trimtags(ised[isat],$
         except=['z','ra','dec']))

; for the centrals need to sort the BCGs file
       bcgs = mrdfits(datapath+'dr8_run_redmapper_v5.2_lgt20_catalog.fits.gz',1)
       match, phot[icen].mem_match_id, bcgs.mem_match_id, m1, m2
       srt = sort(m1)
       bcgs = bcgs[m2[srt]]
;      print, minmax(bcgs.ra-mems[icen].ra), minmax(bcgs.dec-mems[icen].dec)
;      help, where(bcgs.z_lambda-mems[icen].z ne 0)

       centrals = struct_addtags(struct_trimtags(phot[icen],$
         except=['maggies','ivarmaggies']),struct_trimtags(ised[icen],$
         except=['z','ra','dec']))
       centrals = struct_addtags(temporary(centrals),struct_trimtags($
         bcgs,select=['lambda_chisq','lambda_chisq_e','z_lambda','z_lambda_e',$
         'pzbins','pz','bcg_spec_z']))
    endif

return
end
