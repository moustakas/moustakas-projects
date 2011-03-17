function zpt_bootescat2mag, bootes, maggies=maggies, $
  ivarmaggies=ivarmaggies, psf=psf
; simple wrapper to pack the BOOTES photometry in a convenient structure
    bootes_to_maggies, bootes, maggies, ivarmaggies, $
      /nozpoffset, psf=psf
    mag = maggies2mag(maggies,magerr=magerr,ivarmaggies=ivarmaggies)

    result = {$
      ximage:   0.0, $
      yimage:   0.0, $
      bw:       0.0, $
      r:        0.0, $
      i:        0.0, $
      bwerr:    0.0, $
      rerr:     0.0, $
      ierr:     0.0}
    result = replicate(result,n_elements(bootes))
    
    result.ximage = bootes.i_x_image
    result.yimage = bootes.i_y_image

    result.bw = reform(mag[0,*]) & result.bwerr = reform(magerr[0,*])
    result.r = reform(mag[1,*]) & result.rerr = reform(magerr[1,*])
    result.i = reform(mag[2,*]) & result.ierr = reform(magerr[2,*])
    
return, result
end

