pro write_02hicks, outdata
; jm04nov30uofa

    Bmsun = 5.42
    niiha = 0.5D ; Ha = (Ha+[N II]) / (1+[N II]/Ha)
    
    root = '02hicks'
    
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

    data.M_B = raw.M_B - 5.0*alog10(50.0/70.0) ; H_0 = 50 --> 70 km/s/Mpc
    data.M_B_err = 0.1

    lum = 10.0D^(0.4*Bmsun) * 10.0^(-0.4*data.M_B)   ; [L_sun]
    lum_err = alog(10.0) * 0.4D * lum * data.M_B_err ; [L_sun]
       
    data.B_lum_err = lum_err/lum/alog(10.0) ; log L_sun
    data.B_lum = alog10(lum)                ; log L_sun
    
; ---------------------------------------------------------------------------
; fluxes; the H-alpha fluxes must be corrected for contamination from
; [N II]
; ---------------------------------------------------------------------------

    good = where(raw.oii_3727 gt -900,ngood)
    if (ngood ne 0L) then data[good].oii_3727 = 1D-18*transpose([ [raw[good].oii_3727], [raw[good].oii_3727_err] ])

    good = where(raw.h_gamma gt -900,ngood)
    if (ngood ne 0L) then data[good].h_gamma = 1D-18*transpose([ [raw[good].h_gamma], [raw[good].h_gamma_err] ])

    good = where(raw.h_beta gt -900,ngood)
    if (ngood ne 0L) then data[good].h_beta = 1D-18*transpose([ [raw[good].h_beta], [raw[good].h_beta_err] ])

    good = where(raw.oiii_5007 gt -900,ngood)
    if (ngood ne 0L) then begin
       data[good].oiii_5007 = 1D-18*transpose([ [raw[good].oiii_5007], [raw[good].oiii_5007_err] ])
       data[good].oiii_4959 = data[good].oiii_5007 / 3.0
    endif

    good = where(raw.h_alpha gt -900.0,ngood)
    if (ngood ne 0L) then begin
       data[good].h_alpha = 1D-18*transpose([ [raw[good].h_alpha], [raw[good].h_alpha_err] ]) / (1+niiha)
    endif

; ---------------------------------------------------------------------------    
; EW's; assume the error is the same as for the fluxes; the H-alpha
; EW's must be transformed to the rest frame
; ---------------------------------------------------------------------------    

    good = where(raw.oii_3727_ew gt -900,ngood)
    if (ngood ne 0L) then data[good].oii_3727_ew = transpose([ [raw[good].oii_3727_ew], [raw[good].oii_3727_ew*raw[good].oii_3727_err/raw[good].oii_3727] ])

    good = where(raw.h_gamma_ew gt -900,ngood)
    if (ngood ne 0L) then data[good].h_gamma_ew = transpose([ [raw[good].h_gamma_ew], [raw[good].h_gamma_ew*raw[good].h_gamma_err/raw[good].h_gamma] ])

    good = where(raw.h_beta_ew gt -900,ngood)
    if (ngood ne 0L) then data[good].h_beta_ew = transpose([ [raw[good].h_beta_ew], [raw[good].h_beta_ew*raw[good].h_beta_err/raw[good].h_beta] ])

    good = where(raw.oiii_5007_ew gt -900,ngood)
    if (ngood ne 0L) then begin
       data[good].oiii_5007_ew = transpose([ [raw[good].oiii_5007_ew], [raw[good].oiii_5007_ew*raw[good].oiii_5007_err/raw[good].oiii_5007] ])
       data[good].oiii_4959_ew = data[good].oiii_5007_ew / 3.0
    endif

    good = where(raw.h_alpha_ew gt -900,ngood)
    if (ngood ne 0L) then begin
       ew = raw[good].h_alpha_ew / (1+data[good].z_obj)
       ewerr = raw[good].h_alpha_ew*raw[good].h_alpha_err/raw[good].h_alpha / (1+data[good].z_obj)
       data[good].h_alpha_ew = transpose([ [ew], [ewerr] ])
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
