function zpt_zbootescat2mag, zbootes, maggies=maggies, ivarmaggies=ivarmaggies
; simple wrapper to pack the ZBOOTES photometry in a convenient structure
    zbootes_to_maggies, zbootes, maggies, ivarmaggies, /nozpoffset
    mag = maggies2mag(maggies,magerr=magerr,ivarmaggies=ivarmaggies)

    result = {$
      z:        0.0, $
      zerr:     0.0}
;     z_flags:    0}
    result = replicate(result,n_elements(zbootes))

    result.z = reform(mag[0,*]) & result.zerr = reform(magerr[0,*])
;   result.z_flags = zbootes.z_flags

return, result
end

