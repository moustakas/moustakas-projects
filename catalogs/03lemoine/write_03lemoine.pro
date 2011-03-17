pro write_03lemoine, outdata
; jm04dec02uofa

    Bmsun = 5.42

    root = '03lemoine'
    
    path = getenv('CATALOGS_DIR')+'/'+root+'/'
    raw = rsex(path+root+'.dat')
    ngalaxy = n_elements(raw)

    data = init_cat_linefit(ngalaxy=ngalaxy)

; ---------------------------------------------------------------------------    
; properties of interest
; ---------------------------------------------------------------------------    

    data.galaxy = raw.galaxy
    data.z_obj = raw.z
    data.distance = dluminosity(data.z_obj,/Mpc)

    data.M_B = raw.M_B
    data.M_B_err = 0.1

    lum = 10.0D^(0.4*Bmsun) * 10.0^(-0.4*data.M_B)   ; [L_sun]
    lum_err = alog(10.0) * 0.4D * lum * data.M_B_err ; [L_sun]
       
    data.B_lum_err = lum_err/lum/alog(10.0) ; log L_sun
    data.B_lum = alog10(lum)                ; log L_sun
    
; ---------------------------------------------------------------------------
; fluxes    
; ---------------------------------------------------------------------------

    good = where(raw.oii gt -900,ngood)
    if (ngood ne 0L) then data[good].oii_3727 = 1D-17*transpose([ [raw[good].oii], [raw[good].oii_err] ])

    good = where(raw.oiii gt -900,ngood)
    if (ngood ne 0L) then begin
       data[good].oiii_5007 = 1D-17*transpose([ [raw[good].oiii], [raw[good].oiii_err] ])
       data[good].oiii_4959 = data[good].oiii_5007 / 3.0
    endif

    good = where(raw.nii gt -900,ngood)
    if (ngood ne 0L) then begin
       data[good].nii_6584 = 1D-17*transpose([ [raw[good].nii], [raw[good].nii_err] ])
       data[good].nii_6548 = data[good].nii_6584 / 3.0
    endif

    good = where(raw.ha gt -900,ngood)
    if (ngood ne 0L) then data[good].h_alpha = 1D-17*transpose([ [raw[good].ha], [raw[good].ha_err] ])

; ---------------------------------------------------------------------------    
; continuum fluxes
; ---------------------------------------------------------------------------    

; ---------------------------------------------------------------------------    
; compute abundances
; ---------------------------------------------------------------------------    

    adata = im_abundance(data,snrcut=1.0)
    outdata = struct_addtags(data,adata)

    splog, 'Writing '+path+root+'.fits.'
    mwrfits, outdata, path+root+'.fits', /create
    spawn, ['gzip -f ']+path+root+'.fits', /sh

return
end    
