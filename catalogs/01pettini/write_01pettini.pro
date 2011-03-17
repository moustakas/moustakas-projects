pro write_01pettini, outdata
; jm04dec02uofa

    Bmsun = 5.42
    hbeta_abs = 5.0       ; Balmer absorption correction

    root = '01pettini'
    
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

    good = where(raw.oii_3727 gt -900,ngood)
    if (ngood ne 0L) then data[good].oii_3727 = 1D-17*transpose([ [raw[good].oii_3727], [raw[good].oii_3727_err] ])

    good = where(raw.h_beta gt -900,ngood)
    if (ngood ne 0L) then data[good].h_beta = 1D-17*transpose([ [raw[good].h_beta], [raw[good].h_beta_err] ])

    good = where(raw.oiii_5007 gt -900,ngood)
    if (ngood ne 0L) then begin
       data[good].oiii_5007 = 1D-17*transpose([ [raw[good].oiii_5007], [raw[good].oiii_5007_err] ])
       data[good].oiii_4959 = data[good].oiii_5007 / 3.0
    endif

; ---------------------------------------------------------------------------    
; EW's; use the same error as the flux error
; ---------------------------------------------------------------------------    p

    good = where(raw.ew_h_beta gt -900,ngood)
    if (ngood ne 0L) then begin
       data[good].h_beta_ew[0] = raw[good].ew_h_beta
       data[good].h_beta_ew[1] = data[good].h_beta_ew[0] * data[good].h_beta[1] / data[good].h_beta[0]
    endif

; ---------------------------------------------------------------------------    
; continuum fluxes
; ---------------------------------------------------------------------------    

; ---------------------------------------------------------------------------    
; de-redden and compute abundances
; ---------------------------------------------------------------------------    

    adata = im_abundance(data,snrcut=1.0)
    outdata = struct_addtags(data,adata)

    splog, 'Writing '+path+root+'.fits.'
    mwrfits, outdata, path+root+'.fits', /create
    spawn, ['gzip -f ']+path+root+'.fits', /sh

; reddening-corrected fluxes    
    
    idata = iunred_linedust(data,snrcut=1.0,/silent)
    adata = im_abundance(idata,snrcut=1.0)
    outdata = struct_addtags(idata,adata)
    
    splog, 'Writing '+path+root+'_nodust.fits.'
    mwrfits, outdata, path+root+'_nodust.fits', /create
    spawn, ['gzip -f ']+path+root+'_nodust.fits', /sh
    
return
end    
