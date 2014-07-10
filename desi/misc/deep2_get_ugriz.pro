function deep2_get_ugriz, cat, unwise=unwise, degrade=degrade
; jm14may28siena - unpack the multiband photometry of a given a DEEP2
; catalog; also optionally degrade the grz photometry to the
; anticipated depth of our DECam imaging

; compute dust-corrected CFHTLS/ugriz magnitudes
    outcat = cat
    ngal = n_elements(outcat)
    
    deep2_to_maggies, outcat, maggies, ivarmaggies, unwise=unwise, filter=filt

; simplify the input catalog by removing unnecessary information and
; renaming some very annoying structure tags
    outcat = struct_trimtags(outcat,except=['RA_SDSS','DEC_SDSS',$
      'U','G','R','I','Z','UERR','GERR','RERR','IERR','ZERR','W1_*','W2_*'])
    
    oldtags = tag_names(outcat)
    newtags = repstr(repstr(oldtags,'RA_DEEP','ra'),'DEC_DEEP','DEC')
    newtags = repstr(repstr(repstr(newtags,'BESTB','B'),'BESTR','R'),'BESTI','I')
    newtags = repstr(repstr(repstr(newtags,'BESTBERR','BERR'),'BESTRERR','RERR'),$
      'BESTIERR','IERR')
    outcat = im_struct_trimtags(outcat,select=oldtags,newtags=newtags)

; add new tag names    
    outcat = struct_addtags(outcat,replicate({cfhtls_u: -999.0, cfhtls_g: -999.0, $
      cfhtls_r: -999.0, cfhtls_i: -999.0, cfhtls_z: -999.0, cfhtls_uerr: -999.0, $
      cfhtls_gerr: -999.0, cfhtls_rerr: -999.0, cfhtls_ierr: -999.0, cfhtls_zerr: -999.0, $
      w1: -999.0, w2: -999.0, w1err: -999.0, w2err: -999.0},ngal))

; pack the photometry back in; for objects detected at <1-sigma,
; set the uncertainty to zero and the magnitude to the upper limit 
    mag = maggies2mag(maggies,ivarmaggies=ivarmaggies,$
      magerr=magerr,magnsigma=maglimit,nsigma=1.0)

    for ii = 0, n_elements(filt)-1 do begin
       ww = where(mag[ii,*] lt 0 and maglimit[ii,*] gt 0,nww)
       if nww ne 0L then begin
          mag[ii,ww] = maglimit[ii,ww]
          magerr[ii,ww] = 0
       endif
    endfor

    outcat.cfhtls_u = reform(mag[3,*])
    outcat.cfhtls_g = reform(mag[4,*])
    outcat.cfhtls_r = reform(mag[5,*])
    outcat.cfhtls_i = reform(mag[6,*])
    outcat.cfhtls_z = reform(mag[7,*])

    outcat.cfhtls_uerr = reform(magerr[3,*])
    outcat.cfhtls_gerr = reform(magerr[4,*])
    outcat.cfhtls_rerr = reform(magerr[5,*])
    outcat.cfhtls_ierr = reform(magerr[6,*])
    outcat.cfhtls_zerr = reform(magerr[7,*])

    if keyword_set(unwise) then begin
       outcat.w1 = reform(mag[8,*])
       outcat.w2 = reform(mag[9,*])
       outcat.w1err = reform(magerr[8,*])
       outcat.w2err = reform(magerr[9,*])
    endif

; now degrade the grz photometry
    if keyword_set(degrade) then begin
       grzdepth = [24.0, 23.6, 23.0]

       good = where(outcat.cfhtls_gerr gt 0,ngood)
       gflux = 10^((22.5-outcat[good].cfhtls_g)/2.5) + $
         0.16*randomn(1234,ngood)*10^((22.5-grzdepth[0])/2.5)
       outcat[good].cfhtls_g = 22.5 - 2.5*alog10(gflux > 0.01)

       rflux = 10^((22.5-outcat[good].cfhtls_r)/2.5) + $
         0.16*randomn(2345,ngood)*10^((22.5-grzdepth[1])/2.5)
       outcat[good].cfhtls_r = 22.5 - 2.5*alog10(rflux > 0.01)

       zflux = 10^((22.5-outcat[good].cfhtls_z)/2.5) + $
         0.16*randomn(3456,ngood)*10^((22.5-grzdepth[2])/2.5)
       outcat[good].cfhtls_z = 22.5 - 2.5*alog10(zflux > 0.01)
    endif
    
return, outcat
end

