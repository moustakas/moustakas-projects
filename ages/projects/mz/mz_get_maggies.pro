function mz_get_maggies, phot, sdss=sdss, ivarmaggies=ivarmaggies, $
  filterlist=filterlist
; jm10nov01ucsd - simple wrapper to pull out the photometry we want
; for the MZ project; called by BUILD_MZ_PARENT_SAMPLE and MZ_ISEDFIT 

    common com_sdss_photometry, sdssphot, sdsstwomass

    if keyword_set(sdss) then begin
       vagcpath = getenv('VAGC_REDUX')+'/'
       if (n_elements(sdssphot) eq 0L) then sdssphot = $
         mrdfits(vagcpath+'object_sdss_imaging.fits.gz',$
         row=phot.object_position,1)
       if (n_elements(sdsstwomass) eq 0L) then sdsstwomass = $
         mrdfits(vagcpath+'object_twomass.fits.gz',$
         row=phot.object_position,1)

       sdss_to_maggies, smaggies, sivarmaggies, calib=sdssphot
       twomass_to_maggies, sdsstwomass, tmaggies, tivarmaggies
       maggies = [smaggies,tmaggies]
       ivarmaggies = [sivarmaggies,tivarmaggies]
       filterlist = [sdss_filterlist(),twomass_filterlist()]
    endif else begin
       mz_to_maggies, phot, maggies, ivarmaggies, /itot, $
         filterlist=filterlist, use_aper='04', /totalmag
    endelse

return, maggies
end
