pro write_97hammer, outdata
; jm04dec020uofa

    Bmsun = 5.42
    light = 2.99792458D18 ; speed of light [A/s]

    hbeta_abs = 5.0       ; Balmer absorption correction
    halpha_abs = 3.0      ; Balmer absorption correction
    
    root = '97hammer'
    
    path = getenv('CATALOGS_DIR')+'/'+root+'/'
    raw = mrdfits(path+root+'_vizier.fits.gz',1)
    ngalaxy = n_elements(raw)
    
    data = init_cat_linefit(ngalaxy=ngalaxy)

; ---------------------------------------------------------------------------    
; properties of interest
; ---------------------------------------------------------------------------    

    data.galaxy = 'CFRS'+raw.cfrs
    data.z_obj = raw.z
    data.distance = dluminosity(data.z_obj,/Mpc)

    data.M_B = raw.MB + 0.09 - 5.0*alog10(50.0/70.0) ; convert to Vega and H_0 = 50 --> 70 km/s/Mpc
    data.M_B_err = 0.1

    lum = 10.0D^(0.4*Bmsun) * 10.0^(-0.4*data.M_B)   ; [L_sun]
    lum_err = alog(10.0) * 0.4D * lum * data.M_B_err ; [L_sun]
       
    data.B_lum_err = lum_err/lum/alog(10.0) ; log L_sun
    data.B_lum = alog10(lum)                ; log L_sun
    
; ---------------------------------------------------------------------------    
; EW's, which must be transformed into the rest frame; also, H-beta
; must be corrected for stellar absorption; the factor of ten converts
; from the measurement (but not the uncertainty!!!) from units of
; nanometers to Angstroms; consequently, the 9999 and 9998 flags (see
; Hammer et al.) have been divided by 10
; ---------------------------------------------------------------------------    

    good = where((raw.w_oii_ gt 0.0) and (raw.w_oii_ lt 990.0),ngood)
    if (ngood ne 0L) then begin

       data[good].oii_3727_ew[0] = 10.0 * raw[good].w_oii_ / (1+data[good].z_obj)
       data[good].oii_3727_ew[1] = raw[good].e_w_oii_ / (1+data[good].z_obj)

    endif

    good = where((raw.w_oiii_ gt 0.0) and (raw.w_oiii_ lt 990.0),ngood)
    if (ngood ne 0L) then begin

       data[good].oiii_5007_ew[0] = 10.0 * raw[good].w_oiii_ / (1+data[good].z_obj)
       data[good].oiii_5007_ew[1] = raw[good].e_w_oiii_ / (1+data[good].z_obj)

       data[good].oiii_4959_ew = data[good].oiii_5007_ew / 3.0

    endif

    good = where((raw.w_hbeta_ gt 0.0) and (raw.w_hbeta_ lt 990.0),ngood)
    if (ngood ne 0L) then begin

       data[good].h_beta_ew[0] = 10.0 * raw[good].w_hbeta_ / (1+data[good].z_obj) + hbeta_abs
       data[good].h_beta_ew[1] = raw[good].e_w_hbeta_ / (1+data[good].z_obj)

    endif

    good = where((raw.w_halpha_ gt 0.0) and (raw.w_halpha_ lt 990.0),ngood)
    if (ngood ne 0L) then begin

       data[good].h_alpha_ew[0] = 10.0 * raw[good].w_halpha_ / (1+data[good].z_obj) + halpha_abs
       data[good].h_alpha_ew[1] = raw[good].e_w_halpha_ / (1+data[good].z_obj)

    endif

; ---------------------------------------------------------------------------
; fluxes are in units of 10^-29 Angstrom*erg/s/cm2/Hz
; ---------------------------------------------------------------------------

    good = where((raw.oii gt 0.0) and (raw.oii lt 9990.0),ngood)
    if (ngood ne 0L) then begin

       data[good].oii_3727[0] = raw[good].oii*1D-29*light/(data[good].oii_3727_wave*(1+data[good].z_obj))^2
       data[good].oii_3727[1] = raw[good].e_oii*1D-29*light/(data[good].oii_3727_wave*(1+data[good].z_obj))^2

    endif

    good = where((raw.oiii gt 0.0) and (raw.oiii lt 9990.0),ngood)
    if (ngood ne 0L) then begin

       data[good].oiii_5007[0] = raw[good].oiii*1D-29*light/(data[good].oiii_5007_wave*(1+data[good].z_obj))^2
       data[good].oiii_5007[1] = raw[good].e_oiii*1D-29*light/(data[good].oiii_5007_wave*(1+data[good].z_obj))^2

       data[good].oiii_4959 = data[good].oiii_5007 / 3.0

    endif

; to estimate the H-beta flux use the continuum flux at [O III] 5007 
    
    good = where((data.oiii_5007_ew[1] gt 0.0) and (data.oiii_5007[1] gt 0.0) and (data.h_beta_ew[1] gt 0.0),ngood)
    if (ngood ne 0L) then begin

       cflux = data[good].oiii_5007[0] / data[good].oiii_5007_ew[0]
       
       data[good].h_beta[0] = data[good].h_beta_ew[0]*cflux
       data[good].h_beta[1] = data[good].h_beta_ew[1]*cflux

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
