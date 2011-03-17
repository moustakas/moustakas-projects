pro build_oiine_sample, sdss=sdss
; jm09apr01nyu - build the DEEP2 and SDSS OIINE samples 

    oiinepath = deep2_path(/projects)+'oiine/'
    tt = mrdfits(getenv('IMPRO_DIR')+'/im_nebular/temden_table.fits',1)

    snrcut1 = 10.0
    
    if keyword_set(sdss) then begin
; nearly all the relevant sample definitions are determined by the
; choice of the VAGC SAMPLE/LETTER/POSTSTR
       sample = 'dr7'
       letter = 'bsafe'
       poststr = '32'

       ancillary1 = read_sdss_vagc_mpa(/ancillary,sample=sample,$
         letter=letter,poststr=poststr)
       ispec1 = read_sdss_vagc_mpa(/ispec,sample=sample,$
         letter=letter,poststr=poststr)

       sample_zmin = min(ancillary1.z)
       sample_zmax = max(ancillary1.z)
       these = where((ancillary1.z ge sample_zmin) and (ancillary1.z le sample_zmax) and $
         (ispec1.sii_6716[0]/ispec1.sii_6716[1] gt snrcut1) and $
         (ispec1.sii_6731[0]/ispec1.sii_6731[1] gt snrcut1))
       ancillary = ancillary1[these]
       ispec = ispec1[these]

; compute the electron density and errors
       splog, 'Computing electron densities...'
       ratio = ispec.sii_6716[0]/ispec.sii_6731[0] ; use fluxes!!
       ratio_err = im_compute_error(ispec.sii_6716[0],$
         ispec.sii_6716[1],ispec.sii_6731[0],$
         ispec.sii_6731[1],/quotient)

; useful QA checks    
    
;      im_plothist, ispec.sii_6716_sigma[0], bin=5                  
;      im_plothist, ispec[ww].sii_6716_sigma[0], bin=5, /fill, /over
;      im_plothist, ispec.sii_6731_sigma[0], bin=5                  
;      im_plothist, ispec[ww].sii_6731_sigma[0], bin=5, /fill, /over
;   
;      im_plothist, ispec.sii_6716_chi2, bin=0.1, xr=[0,10]                 
;      im_plothist, ispec[ww].sii_6716_chi2, bin=0.1, /fill, /over
;      im_plothist, ispec.sii_6731_chi2, bin=0.1, xr=[0,10]
;      im_plothist, ispec[ww].sii_6731_chi2, bin=0.1, /fill, /over

       tt0 = systime(1)
       temden = im_temden('s_ii',ratio,ratio_err=ratio_err,$
         nmonte=0L,temp_guess=1E4,dens_guess=1E2)
       splog, 'Total time = '+strtrim(string((systime(1)-$
         tt0)/60.0,format='(F12.1)'),2)+' minutes.'
    
; write out       
       out = struct_addtags(temden,ispec)
       im_mwrfits, out, oiinepath+'sdss_oiine_ispec.fits'
       im_mwrfits, ancillary, oiinepath+'sdss_oiine_ancillary.fits'

    endif else begin
       ancillary1 = read_deep2(/ancillary)
       ispec1 = read_deep2(/ispec)

       sample_zmin = 0.7
       sample_zmax = 1.5
       these = where((ancillary1.z ge sample_zmin) and (ancillary1.z le sample_zmax) and $
         (ispec1.oii_3727_1[0]/ispec1.oii_3727_1[1] gt snrcut1) and $
         (ispec1.oii_3727_2[0]/ispec1.oii_3727_2[1] gt snrcut1))
       ancillary = ancillary1[these]
       ispec = ispec1[these]

; compute the electron density and errors
       splog, 'Computing electron densities...'
       ratio = 1.02*ispec.oii_3727_1_ew[0]/ispec.oii_3727_2_ew[0] ; corrected ratio
       ratio_err = im_compute_error(ispec.oii_3727_1_ew[0],$
         ispec.oii_3727_1_ew[1],ispec.oii_3727_2_ew[0],$
         ispec.oii_3727_2_ew[1],/quotient)

; useful QA checks    
    
;      ww = where(ratio lt 0.67) & help, ww
;      ww = where(ratio+ratio_err lt 0.67) & help, ww
;   
;      im_plothist, ispec.oii_3727_1_sigma[0], bin=5                  
;      im_plothist, ispec[ww].oii_3727_1_sigma[0], bin=5, /fill, /over
;      im_plothist, ispec.oii_3727_2_sigma[0], bin=5                  
;      im_plothist, ispec[ww].oii_3727_2_sigma[0], bin=5, /fill, /over
;   
;      im_plothist, ispec.oii_3727_1_chi2, bin=0.1, xr=[0,10]                 
;      im_plothist, ispec[ww].oii_3727_1_chi2, bin=0.1, /fill, /over
;      im_plothist, ispec.oii_3727_2_chi2, bin=0.1, xr=[0,10]
;      im_plothist, ispec[ww].oii_3727_2_chi2, bin=0.1, /fill, /over

       tt0 = systime(1)
       temden = im_temden('o_ii_dens',ratio,ratio_err=ratio_err,$
         nmonte=0L,temp_guess=1E4,dens_guess=1E2)
       splog, 'Total time = '+strtrim(string((systime(1)-$
         tt0)/60.0,format='(F12.1)'),2)+' minutes'

; write out       
       out = struct_addtags(temden,ispec)
       im_mwrfits, out, oiinepath+'deep_oiine_ispec.fits'
       im_mwrfits, ancillary, oiinepath+'deep_oiine_ancillary.fits'

    endelse    

return
end
    
