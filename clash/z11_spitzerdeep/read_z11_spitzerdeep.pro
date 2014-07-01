function read_z11_spitzerdeep
; jm14jun23siena - based on READ_Z11() but with updated Spitzer
; photometry 

    path = getenv('CLASH_PROJECTS')+'/z11_spitzerdeep/'
    cat = rsex(path+'M0647JD_final2sum.cat')

; pick out JD2 (from Dan: "Yes, and the 2nd row of data (id 20000) is
; our primary target JD2, the 2nd-brightest but still fairly well
; isolated image.")
    ww = where(cat.id eq 20000)
    cat = cat[ww]

; add the redshifts, magnifications, and updated Spitzer photometry 
    cat = struct_addtags(cat,{z: 10.6, mu: 1.0, galaxy: ''})
    cat.mu = 7.0
    cat.galaxy = 'JD2'

; replicate so we can fit the two groups that are doing IRAC
; photometry; all fluxes are in microJy 
    cat = replicate(cat,3)
    
; Egami's photometry     
    cat[0].galaxy = 'JD2 (Egami)'
    cat[0].ch1_flux = 0
    cat[0].ch1_fluxerr = 166*1E-3
    cat[0].ch2_flux = 436*1E-3
    cat[0].ch2_fluxerr = 138*1E-3
    
; Xinwen's photometry    
    cat[1].galaxy = 'JD2 (Xinwen)'
    cat[1].ch1_flux = 0.0*1E-3    
    cat[1].ch1_fluxerr = 114*1E-3 
    cat[1].ch2_flux = 158*1E-3    
    cat[1].ch2_fluxerr = 71*1E-3  

; Renske's photometry     
    cat[2].galaxy = 'JD2 (Renske)'
    cat[2].ch1_flux = 291.9*1E-3
    cat[2].ch1_fluxerr = 226.1*1E-3
    cat[2].ch2_flux = 185.1*1E-3
    cat[2].ch2_fluxerr = 107.4*1E-3
    
return, cat
end    
