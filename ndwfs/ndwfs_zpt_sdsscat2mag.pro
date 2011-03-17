function ndwfs_zpt_sdsscat2mag, sdss, flux=flux, maggies=maggies, $
  ivarmaggies=ivarmaggies
; simple wrapper to pack the SDSS photometry in a convenient structure
    if (n_elements(flux) eq 0) then flux = 'psf'
    sdss_to_maggies, maggies, ivarmaggies, calib=sdss, flux=flux
    mag = maggies2mag(maggies,ivarmaggies=ivarmaggies,magerr=magerr)
    
    result = {$
      u: 0.0, $
      g: 0.0, $
      r: 0.0, $
      i: 0.0, $
      z: 0.0, $
      uerr: 0.0, $
      gerr: 0.0, $
      rerr: 0.0, $
      ierr: 0.0, $
      zerr: 0.0}
    result = replicate(result,n_elements(sdss))
    
    result.u = reform(mag[0,*]) & result.uerr = reform(magerr[0,*])
    result.g = reform(mag[1,*]) & result.gerr = reform(magerr[1,*])
    result.r = reform(mag[2,*]) & result.rerr = reform(magerr[2,*])
    result.i = reform(mag[3,*]) & result.ierr = reform(magerr[3,*])
    result.z = reform(mag[4,*]) & result.zerr = reform(magerr[4,*])
    
return, result
end

