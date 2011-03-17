pro write_99glazebrook, outdata
; jm04dec31uofa

    Bmsun = 5.42

; read and cross-match against the Hammer et al. (1997) catalog 

    root = '99glazebrook'
    
    path = getenv('CATALOGS_DIR')+'/'+root+'/'
    raw = rsex(path+root+'.dat')
    ngalaxy = n_elements(raw)

    data = init_cat_linefit(ngalaxy=ngalaxy)

    hammer = read_97hammer()

; ---------------------------------------------------------------------------    
; properties of interest
; ---------------------------------------------------------------------------    

    data.galaxy = raw.cfrsid
    data.z_obj = raw.z
    data.distance = dluminosity(data.z_obj,/Mpc)

;   red, h100=0.5, omega_lambda=0.0, omega_0=1.0
;   g99_distance = dluminosity(data.z_obj,/Mpc)
;   red, h100=0.7, omega_lambda=0.7, omega_0=0.3

    match, data.galaxy, hammer.galaxy, datamatch, hammermatch
    nmatch = n_elements(datamatch)
    niceprint, data[datamatch].galaxy, hammer[hammermatch].galaxy, raw[datamatch].oii_3727_ew, $
      hammer[hammermatch].oii_3727_ew[0], hammer[hammermatch].oii_3727[0]*1D17
;   help, hammermatch

    bigindx = lindgen(ngalaxy)
    remove, datamatch, bigindx
    splog, 'Missing objects:'
    niceprint, data[bigindx].galaxy ; these objects are not in Hammer et al. (1997)
    
    data[datamatch].M_B = hammer[hammermatch].M_B
    data[datamatch].M_B_err = hammer[hammermatch].M_B_err

;   data.M_B = raw.M_B_AB + 0.09 - 5.0*alog10(50.0/70.0) ; convert to Vega and H_0 = 50 --> 70 km/s/Mpc
;   data.M_B_err = 0.1

    good = where(data.m_b gt -900.0,ngood)
    if (ngood ne 0L) then begin
       lum = 10.0D^(0.4*Bmsun) * 10.0^(-0.4*data[good].M_B) ; [L_sun]
       lum_err = alog(10.0) * 0.4D * lum * data[good].M_B_err ; [L_sun]       
       data[good].B_lum_err = lum_err/lum/alog(10.0) ; log L_sun
       data[good].B_lum = alog10(lum) ; log L_sun
    endif
    
; ---------------------------------------------------------------------------
; fluxes    
; ---------------------------------------------------------------------------

    good = where(raw.h_alpha gt -900,ngood)
    if (ngood ne 0L) then begin

       data[good].h_alpha = 1D-17*transpose([ [raw[good].h_alpha], [raw[good].h_alpha_err] ])

       factor = raw[good].ha_lum*1D41/data[good].h_alpha[0]
       
       data[good].oii_3727[0] = raw[good].oii_lum*1D41/factor
       data[good].oii_3727[1] = raw[good].oii_lum_err*1D41/factor

;      niceprint, data[good].h_alpha[0]/data[good].oii_3727[0], raw[good].ha_lum/raw[good].oii_lum

    endif

;   factor = 4.0*!dpi*(g99_distance*1D6*3.086D18)^2
;   oii = raw.oii_lum*1D41/factor
;   oii_err = raw.oii_lum_err*1D41/factor
;   data.oii_3727 = transpose([ [oii], [oii_err] ])
    
; check our luminosity business    
    
;   niceprint, raw.ha_lum*1D41/factor, data.h_alpha[0], raw.ha_lum_err*1D41/factor, data.h_alpha[1]

; don't use the Hammer fluxes! (see comments in 99glazebrook.dat)
    
;   data[datamatch].oii_3727 = hammer[hammermatch].oii_3727
;   data[datamatch].h_beta = hammer[hammermatch].h_beta
;   data[datamatch].oiii_4959 = hammer[hammermatch].oiii_4959
;   data[datamatch].oiii_5007 = hammer[hammermatch].oiii_5007
    
; ---------------------------------------------------------------------------    
; EW's
; ---------------------------------------------------------------------------    

; since the tabulated [O II] EWs are identical use the Hammer data 
    
;   good = where(raw.oii_3727_ew gt -900,ngood)
;   if (ngood ne 0L) then data[good].oii_3727_ew = transpose([ [raw[good].oii_3727_ew], $
;     [raw[good].oii_3727_ew*raw[good].oii_3727_err/raw[good].oii_3727] ])

;   good = where(raw.h_alpha_ew gt -900,ngood)
;   if (ngood ne 0L) then data[good].h_alpha_ew = transpose([ [raw[good].h_alpha_ew], $
;     [raw[good].h_alpha_ew*raw[good].h_alpha_err/raw[good].h_alpha] ])

    data[datamatch].oii_3727_ew = hammer[hammermatch].oii_3727_ew
;   data[datamatch].h_beta_ew = hammer[hammermatch].h_beta_ew
;   data[datamatch].oiii_4959_ew = hammer[hammermatch].oiii_4959_ew
;   data[datamatch].oiii_5007_ew = hammer[hammermatch].oiii_5007_ew

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
