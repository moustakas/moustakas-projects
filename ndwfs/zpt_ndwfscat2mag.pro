function zpt_ndwfscat2mag, ndwfs, maggies=maggies, ivarmaggies=ivarmaggies
; simple wrapper to pack the NDWFS photometry in a convenient structure
    ndwfs_to_maggies, ndwfs, maggies, ivarmaggies, /nozpoffset
    mag = maggies2mag(maggies,magerr=magerr,ivarmaggies=ivarmaggies)

    result = {$
      bw:       0.0, $
      r:        0.0, $
      i:        0.0, $
      k:        0.0, $
      bwerr:    0.0, $
      rerr:     0.0, $
      ierr:     0.0, $
      kerr:     0.0, $
      bw_flags:   0, $
      r_flags:    0, $
      i_flags:    0, $
      k_flags:    0}
    result = replicate(result,n_elements(ndwfs))
    
    result.bw = reform(mag[0,*]) & result.bwerr = reform(magerr[0,*])
    result.r = reform(mag[1,*]) & result.rerr = reform(magerr[1,*])
    result.i = reform(mag[2,*]) & result.ierr = reform(magerr[2,*])
    result.k = reform(mag[3,*]) & result.kerr = reform(magerr[3,*])
    
    result.bw_flags = ndwfs.bw_flags
    result.r_flags = ndwfs.r_flags
    result.i_flags = ndwfs.i_flags
    result.k_flags = ndwfs.k_flags

return, result
end

