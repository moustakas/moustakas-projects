pro write_02tresse, outdata
; jm04dec02uofa

; read and cross-match against the Hammer et al. (1997) catalog 

    Bmsun = 5.42

    root = '02tresse'
    
    path = getenv('CATALOGS_DIR')+'/'+root+'/'
    raw = rsex(path+root+'.dat')
    ngalaxy = n_elements(raw)

    hammer = read_97hammer()

    data = init_cat_linefit(ngalaxy=ngalaxy)

; ---------------------------------------------------------------------------    
; properties of interest
; ---------------------------------------------------------------------------    

    data.galaxy = raw.cfrsid
    data.z_obj = raw.z
    data.distance = dluminosity(data.z_obj,/Mpc)

    data.M_B = raw.M_B_AB + 0.09 - 5.0*alog10(50.0/70.0) ; convert to Vega and H_0 = 50 --> 70 km/s/Mpc
    data.M_B_err = 0.1

    lum = 10.0D^(0.4*Bmsun) * 10.0^(-0.4*data.M_B)   ; [L_sun]
    lum_err = alog(10.0) * 0.4D * lum * data.M_B_err ; [L_sun]
       
    data.B_lum_err = lum_err/lum/alog(10.0) ; log L_sun
    data.B_lum = alog10(lum)                ; log L_sun
    
    match, data.galaxy, hammer.galaxy, datamatch, hammermatch
    niceprint, data[datamatch].galaxy, hammer[hammermatch].galaxy, data[datamatch].M_B, hammer[hammermatch].M_B
    help, hammermatch

    bigindx = lindgen(ngalaxy)
    remove, datamatch, bigindx
    niceprint, data[bigindx].galaxy
    
; ---------------------------------------------------------------------------
; fluxes    
; ---------------------------------------------------------------------------

    good = where(raw.oii_3727 gt -900,ngood)
    if (ngood ne 0L) then data[good].oii_3727 = 1D-17*transpose([ [raw[good].oii_3727], [raw[good].oii_3727_err] ])

    good = where(raw.h_alpha gt -900,ngood)
    if (ngood ne 0L) then data[good].h_alpha = 1D-17*transpose([ [raw[good].h_alpha], [raw[good].h_alpha_err] ])

    data[datamatch].h_beta = hammer[hammermatch].h_beta
    data[datamatch].oiii_4959 = hammer[hammermatch].oiii_4959
    data[datamatch].oiii_5007 = hammer[hammermatch].oiii_5007
    
; ---------------------------------------------------------------------------    
; EW's; assume the error is the same as for the fluxes
; ---------------------------------------------------------------------------    

    good = where(raw.oii_3727_ew gt -900,ngood)
    if (ngood ne 0L) then data[good].oii_3727_ew = transpose([ [raw[good].oii_3727_ew], $
      [raw[good].oii_3727_ew*raw[good].oii_3727_err/raw[good].oii_3727] ])

    good = where(raw.h_alpha_ew gt -900,ngood)
    if (ngood ne 0L) then data[good].h_alpha_ew = transpose([ [raw[good].h_alpha_ew], $
      [raw[good].h_alpha_ew*raw[good].h_alpha_err/raw[good].h_alpha] ])

    data[datamatch].h_beta_ew = hammer[hammermatch].h_beta_ew
    data[datamatch].oiii_4959_ew = hammer[hammermatch].oiii_4959_ew
    data[datamatch].oiii_5007_ew = hammer[hammermatch].oiii_5007_ew

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
