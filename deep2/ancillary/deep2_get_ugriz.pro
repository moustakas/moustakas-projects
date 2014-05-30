function deep2_get_ugriz, cat, unwise=unwise
; jm14may28siena - compute dust-corrected ugriz magnitudes
    deep2_to_maggies, cat, mm, ii, unwise=unwise, filter=filt
    cat = struct_addtags(cat,replicate({ugriz: fltarr(5), $
      ugriz_err: fltarr(5), bri: fltarr(3), bri_err: fltarr(3), $
      wise: fltarr(2), wise_err: fltarr(2)},$
      n_elements(cat)))
    cat = im_struct_trimtags(cat,select=tag_names(cat),newtags=$
      repstr(repstr(tag_names(cat),'RA_DEEP','ra'),'DEC_DEEP','DEC'))
    cat = struct_trimtags(cat,except=['RA_SDSS','DEC_SDSS'])
    mag = maggies2mag(mm,ivarmaggies=ii,magerr=magerr,magnsigma=maglimit)

    for ii = 0, n_elements(filt)-1 do begin
       ww = where(mag[ii,*] lt 0 and maglimit[ii,*] gt 0,nww)
       if nww ne 0L then begin
          mag[ii,ww] = maglimit[ii,ww]
          magerr[ii,ww] = 0
       endif
    endfor
    
    cat.bri = mag[0:2,*]
    cat.bri_err = magerr[0:2,*]
    cat.ugriz = mag[3:7,*]
    cat.ugriz_err = magerr[3:7,*]
    if keyword_set(unwise) then begin
       cat.wise = mag[8:9,*]
       cat.wise_err = magerr[8:9,*]
    endif
return, cat
end

