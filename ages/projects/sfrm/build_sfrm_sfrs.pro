function init_sfrs_table, ngal
; initialize the data table
    
    sfrs = {$
      mips:           0,$ ; is this a MIPS source?
      lir_ce01:     0.0,$ ; Chary & Elbaz (2001) L(IR)
      lir_ce01_err: 0.0,$
      lir_dh02:     0.0,$ ; Dale & Helou (2002) L(IR)
      lir_dh02_err: 0.0,$
      lir_r09:      0.0,$ ; Rieke et al. (2009) L(IR)
      lir_r09_err:  0.0,$

      l24_ce01:     0.0,$ ; Chary & Elbaz (2001) L(24)
      l24_ce01_err: 0.0,$
      l24_dh02:     0.0,$ ; Dale & Helou (2002) L(24)
      l24_dh02_err: 0.0,$
      l24_r09:      0.0,$ ; Rieke et al. (2009) L(24)
      l24_r09_err:  0.0}
    sfrs = replicate(sfrs,ngal)

return, sfrs
end

pro build_sfrm_sfrs, sfrs, sdss=sdss
; jm10feb17ucsd - build the SFR estimates for the AGES/SFRM samples 

    sfrmpath = ages_path(/projects)+'sfrm/'

    if keyword_set(sdss) then begin
       splog, '#########################'
       splog, 'Building the SDSS SFRs'
       sample = read_sfrm_sample(/sdss)
       ngal = n_elements(sample)

; write out       
;      im_mwrfits, sample, sfrmpath+'sdss_sfrs.fits', /clobber
       
    endif else begin
       splog, '#########################'
       splog, 'Building the AGES SFRs'
       sample = read_sfrm_sample()
       ngal = n_elements(sample)

       sfrs = init_sfrs_table(ngal)

; --------------------------------------------------
; compute L(IR) using several different methods and the observed
; 24-micron flux densities; assign undetected objects the 80%
; completeness limit, 0.27mJy
       bright = where(sample.phot_mips24 gt 0.0,comp=faint)
       sfrs[bright].mips = 1
       f24 = sample.phot_mips24
       f24_err = sample.phot_mips24_err
       f24[faint] = 0.27 ; [mJy]

; Chary & Elbaz (2001)       
       sfrs.lir_ce01 = im_f24tolir(sample.z,f24,f24_err=f24_err,/chary,$
         err_lir=lir_ce01_err,l24=l24_ce01,err_l24=l24_ce01_err)
       sfrs.l24_ce01 = l24_ce01
       sfrs.lir_ce01_err = lir_ce01_err
       sfrs.l24_ce01_err = l24_ce01_err
; Dale & Helou (2002)
       sfrs.lir_dh02 = im_f24tolir(sample.z,f24,f24_err=f24_err,/dale,$
         err_lir=lir_dh02_err,l24=l24_dh02,err_l24=l24_dh02_err)
       sfrs.l24_dh02 = l24_dh02
       sfrs.lir_dh02_err = lir_dh02_err
       sfrs.l24_dh02_err = l24_dh02_err
; Rieke et al. (2009)
       sfrs.lir_r09 = im_f24tolir(sample.z,f24,f24_err=f24_err,/rieke,$
         err_lir=lir_r09_err,l24=l24_r09,err_l24=l24_r09_err)
       sfrs.l24_r09 = l24_r09
       sfrs.lir_r09_err = lir_r09_err
       sfrs.l24_r09_err = l24_r09_err

       
       
       
; push these QAplots to a separate routine       
       plot, alog10(sfrs.lir_ce01), alog10(sfrs.lir_ce01/sfrs.lir_dh02), psym=3, xsty=3, ysty=3, yr=[-0.5,0.5]
       plot, alog10(sfrs.lir_ce01), alog10(sfrs.lir_ce01/sfrs.lir_r09), psym=3, xsty=3, ysty=3, yr=[-0.5,0.5]
       plot, sample.z, alog10(sfrs.lir_ce01/sfrs.lir_dh02), psym=3, xsty=3, ysty=3, yr=[-0.5,0.5]
       plot, sample.z, alog10(sfrs.lir_ce01/sfrs.lir_r09), psym=3, xsty=3, ysty=3, yr=[-0.5,0.5]

       plot, alog10(sfrs.l24_ce01), alog10(sfrs.l24_ce01/sfrs.l24_dh02), psym=3, xsty=3, ysty=3, yr=[-0.5,0.5]
       plot, alog10(sfrs.l24_ce01), alog10(sfrs.l24_ce01/sfrs.l24_r09), psym=3, xsty=3, ysty=3, yr=[-0.5,0.5]
       plot, sample.z, alog10(sfrs.l24_ce01/sfrs.l24_dh02), psym=3, xsty=3, ysty=3, yr=[-0.5,0.5]
       plot, sample.z, alog10(sfrs.l24_ce01/sfrs.l24_r09), psym=3, xsty=3, ysty=3, yr=[-0.5,0.5]

stop       

; write out       
       im_mwrfits, sfrs, sfrmpath+'ages_sfrs.fits', /clobber

    endelse
       
return
end
