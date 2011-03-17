pro write_05maier, outdata
; jm05jan01uofa - these objects are a subset of Lilly et al. (2003),
;                 so we merge them there; also, the EWs are more
;                 complete in write_03lilly

    snrcut = 3.0
    Bmsun = k_solar_magnitudes(filterlist='bessell_B.par')
    Bvega2ab = (k_vega2ab(filterlist='bessell_B.par'))[0]
    
    root = '05maier'
    path = getenv('CATALOGS_DIR')+'/'+root+'/'
    raw = rsex(path+root+'.dat')
    ngalaxy = n_elements(raw)

; initialize the output data structure    

    data = init_cat_linefit(ngalaxy=ngalaxy)
    
; ---------------------------------------------------------------------------    
; store relevant galaxy properties
; ---------------------------------------------------------------------------    

    data.id = string(lindgen(ngalaxy),format='(I2.2)')
    data.galaxy = strtrim(raw.galaxy,2)
    data.z_obj = raw.redshift
    data.distance = dluminosity(data.z_obj,/Mpc)

    data.M_B = raw.M_B_AB - Bvega2ab ; convert to Vega
    data.M_B_err = 0.1 ; reference?

    lum = 10.0D^(0.4*Bmsun) * 10.0^(-0.4*data.M_B)   ; [L_sun]
    lum_err = alog(10.0) * 0.4D * lum * data.M_B_err ; [L_sun]
       
    data.B_lum_err = lum_err/lum/alog(10.0) ; log L_sun
    data.B_lum = alog10(lum)                ; log L_sun

; ---------------------------------------------------------------------------
; fluxes    
; ---------------------------------------------------------------------------
    
    good = where(raw.oii gt 0.0)
    data[good].oii_3727[0] = raw[good].oii*1D-17
    data[good].oii_3727[1] = raw[good].oii_err*1D-17

    good = where(raw.hb gt 0.0)
    data[good].h_beta[0] = raw[good].hb*1D-17 ; absorption-corrected
    data[good].h_beta[1] = raw[good].hb_err*1D-17

    good = where(raw.oiii gt 0.0)
    data[good].oiii_5007[0] = raw[good].oiii*1D-17
    data[good].oiii_5007[1] = raw[good].oiii_err*1D-17
    data[good].oiii_4959 = data[good].oiii_5007 / 3.0 ; assume this is true

    good = where(raw.ha gt 0.0)
    data[good].h_alpha[0] = raw[good].ha*1D-17
    data[good].h_alpha[1] = raw[good].ha_err*1D-17

    good = where(raw.nii gt 0.0)
    data[good].nii_6584[0] = raw[good].nii*1D-17
    data[good].nii_6584[1] = raw[good].nii_err*1D-17
    data[good].nii_6548 = data[good].nii_6584 / 3.0 ; assume this is true

; ---------------------------------------------------------------------------    
; compute miscellaneous properties and write out 
; ---------------------------------------------------------------------------    

    adata = im_abundance(data,snrcut=snrcut)
    outdata = struct_addtags(data,adata)
    outdata = abundance_catalogs_log12oh(outdata)

    lums = compute_linelums(outdata,select_lines=['oii_3727','oiii_5007','h_beta','h_alpha'])
    outdata = struct_addtags(outdata,lums)
    sfrs = compute_sfrs(outdata)

    class = iclassification(outdata,/kauffmann)
    outdata = struct_addtags(struct_addtags(outdata,sfrs),class)

    splog, 'Writing '+path+root+'.fits.'
    mwrfits, outdata, path+root+'.fits', /create
    spawn, ['gzip -f ']+path+root+'.fits', /sh

; ---------------------------------------------------------------------------    
; correct for reddening, compute miscellaneous properties and write out 
; ---------------------------------------------------------------------------    

    idata = iunred_linedust(data,snrcut=snrcut,/silent,/nopropagate)
    adata = im_abundance(idata,snrcut=snrcut)
    outdata = struct_addtags(idata,adata)
    outdata = abundance_catalogs_log12oh(outdata)

    lums = compute_linelums(outdata,select_lines=['oii_3727','oiii_5007','h_beta','h_alpha'])
    outdata = struct_addtags(outdata,lums)
    sfrs = compute_sfrs(outdata)

    class = iclassification(outdata,/kauffmann)
    outdata = struct_addtags(struct_addtags(outdata,sfrs),class)

    splog, 'Writing '+path+root+'_nodust.fits.'
    mwrfits, outdata, path+root+'_nodust.fits', /create
    spawn, ['gzip -f ']+path+root+'_nodust.fits', /sh
    
return
end    
