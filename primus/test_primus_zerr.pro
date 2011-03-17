pro test_primus_zerr, rerun, mask=mask
; jm09jul21nyu 

; for this rerun, figure out which masks have been reduced
    basedir = getenv('PRIMUS_DATA')+'/redux/1d/'+rerun+'/'
    zallfiles = file_search(basedir+'ut??????/*-zAll.fits.gz')
    allmasks = strmid(file_basename(zallfiles),0,8)
    if (n_elements(mask) ne 0) then begin
       this = where(allmasks eq mask,nthis)
       if (nthis eq 0) then begin
          splog, 'Mask '+mask+' not found'
          return
       endif
       allmasks = allmasks[this]
       zallfiles = zallfiles[this]
    endif
    nmask = n_elements(allmasks)

    for imask = 0L, nmask-1L do begin

       splog, 'Reading '+file_basename(zallfiles[imask])
       oned1 = mrdfits(zallfiles[imask],1)

       good = where(oned1.knownz gt 0.0,ngood)
       if (ngood ne 0L) then begin
          oned = oned1[good]
          
       endif

stop    

; this code is crap/testing:    
;    
;; interpolate the PDF such that ZMIN is centered on a pixel    
;    zmindex = primus_findex(osamp_z,zmin*1.0D) ; [pixel]
;    zanchor = (round(zmindex)<(osamp_npix-1L))>0L
;    new_index = osamp_index+zmindex-zanchor ; recenter on ZMIN
;    quadterp, osamp_index, osamp_pdf, new_index, new_pdf, missing=0.0
;    quadterp, osamp_index, osamp_z, new_index, new_z, missing=0.0
;;   linterp, osamp_index, osamp_pdf, new_index, new_pdf, missing=0.0
;;   linterp, osamp_index, osamp_z, new_index, new_z, missing=0.0
;    new_pdf = new_pdf/total(new_pdf)
;
;    if keyword_set(debug) then begin
;       ppdf = exp(-0.5*(chi2-chi2min))
;       ppdf = ppdf/total(ppdf)
;       djs_plot, z, ppdf, xsty=3, ysty=3, psym=-4, xr=zmin+[-0.02,+0.02]
;       djs_oplot, osamp_z, osamp_pdf, color='orange', psym=-4
;       djs_oplot, new_z, new_pdf, color='red', psym=-4
;       djs_oplot, zmin*[1,1], !y.crange, line=0
;    endif
;
;    print, interpol(new_index,new_z,zmin)
;    
;; intergrate the PDF to higher redshift
;    newpdfhi = newpdf
;    newpdfhi[0L:zanchor-1L] = 0.0
;    newpdfhi = newpdfhi/total(newpdfhi)
;    zhi = interpol(z,total(newpdfhi,/cumu),sigma)
;
;; ...and to lower redshift       
;    newpdflo = newpdf
;    newpdflo[zanchor+1L:npix-1L] = 0.0
;    newpdflo = newpdflo/total(newpdflo)
;    zlo = interpol(z,total(newpdflo,/cumu),sigma)
;    
;    pixsigma = weighted_quantile(newindex,newpdf,$
;      quant=[0.5-sigma/2.0D,sigma/2.0+0.5])
;    zsigma = interpol(osamp_z,osamp_index,pixsigma)
;    zhierr = (zsigma[1]-zmin)/3.0 ; 1-sigma
;    zloerr = (zmin-zsigma[0])/3.0 ; 1-sigma
;    zerr = djs_mean([zloerr,zhierr])


return
end   
