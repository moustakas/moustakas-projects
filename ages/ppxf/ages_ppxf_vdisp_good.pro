function ages_ppxf_vdisp_good, zabsvdisp, debug=debug
; jm09dec08ucsd - not all the VDISP values computed by PPXF are
;   created equal: the quality of the results depends on a lot of
;   things such as the spectrum S/N and the rest-frame spectral range
;   of the data; this routine, therefore, decides which VDISP
;   measurements are reliable and returns a fiducial VDISP value for
;   crap measurements

; ZABSVDISP is the output from AGES_GET_ZABS_VDISP

    vdisp_fiducial = 165.0 ; [km/s]
    snrcut = 2.0

    vv = zabsvdisp.vdisp
    vverr = zabsvdisp.vdisp_err
    vvsnr = vv/(vverr+(vverr eq 0))*(vverr ne 0)
    vdisp = vv

    crap = where(vvsnr lt snrcut,ncrap,$
      comp=good,ncomp=ngood)
;   if (ngood ne 0) then splog, median(vdisp[good])
    if (ncrap ne 0) then vdisp[crap] = vdisp_fiducial

    if keyword_set(debug) then begin
       djs_plot, zabsvdisp.snr, vvsnr>1E-2, psym=3, /xlog, /ylog, $
         xsty=3, ysty=3, yr=[0.01,100], xtitle='Continuum S/N (pixel^{-1})', $
         ytitle='\sigma/\delta\sigma'
       if (ncrap ne 0) then djs_oplot, zabsvdisp[crap].snr, $
         vvsnr[crap]>1E-2, psym=3, color='red'
       cc = get_kbrd(1)
       djs_plot, zabsvdisp.snr, vv, psym=3, /xlog, /ylog, $
         xsty=3, ysty=3, xtitle='Continuum S/N (pixel^{-1})', $
         ytitle='\sigma (km s^{-1})'
       if (ncrap ne 0) then djs_oplot, zabsvdisp[crap].snr, $
         vv[crap], psym=3, color='red'
    endif

return, vdisp
end    
