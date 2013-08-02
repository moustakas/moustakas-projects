function read_z11
; jm13jul30siena

; from D. Coe, in preparation for our Spitzer proposal:

; As a first test, you might fix z = 10.8 and try the best fitting SED
; values, roughly: 
; ch1 ~ 200 nJy
; ch2 ~ 300 nJy
; 
; with uncertainties that shrink to:
; ch1 ~ 44 nJy
; ch2 ~ 52 nJy
; in 50-hour exposures
; 
; compared to the currently measured uncertainties in 5 hours:
; ch1 ~ 167 nJy
; ch2 ~ 139 nJy
; (larger than expected from SENS-PET, but probably more realistic I
; think...) 
; 
; I wonder how well a Balmer break could be disentangled from [OII]
; emission contributing to the ch2 flux? 
; 
; Other interesting cases might be fluxes of ~170 nJy in both
; channels (~ flat in F_nu) or ~100 nJy in both channels (young
; and pristine, as expected?).  But I'm happy to focus on the
; case above... 
; 
; Interested to hear what you find! 
; 
; Thanks again,
; Dan

    path = getenv('IM_PROJECTS_DIR')+'/clash/z11/'
    cat = rsex(path+'M0647JD_final2sum.cat')

; add the redshift
    cat = struct_addtags(cat,replicate({z: 10.8, logmu: 0.0},n_elements(cat)))

; pick out JD2 (from Dan: "Yes, and the 2nd row of data (id 20000) is
; our primary target JD2, the 2nd-brightest but still fairly well
; isolated image.")
    ww = where(cat.id eq 20000)
    cat = cat[ww]
;   cat.ch1_flux = 0.0
;   cat.ch1_flux = 250.0*1E-3   ; uJy
;   cat.ch2_flux = 436.0*1E-3   ; uJy
;   cat.ch1_flux = 200.0*1E-3 ; uJy
;   cat.ch2_flux = 300.0*1E-3 ; uJy
    
; try two cases: 5 hours of IRAC vs 50 hours
    cat = replicate(cat,2)
    cat.ch1_flux = [0.0,200.0]*1E-3   ; uJy
    cat.ch2_flux = [436.0,350.0]*1E-3   ; uJy
    cat.ch1_fluxerr = [166.0,20.4]*1E-3 ; uJy
    cat.ch2_fluxerr = [139.0,29.4]*1E-3 ; uJy

;; try three cases: 5 hours of IRAC vs 50 hours vs 200 hours
;    cat = replicate(cat,3)
;    cat.ch1_fluxerr = [167.0,20.4,20.4/sqrt(4)]*1E-3 ; uJy
;    cat.ch2_fluxerr = [139.0,29.4,29.4/sqrt(4)]*1E-3 ; uJy

return, cat
end    
