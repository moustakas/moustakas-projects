pro write_92kennicutt, data
; jm04dec31uofa

    Bmsun = 5.42
    light = 2.99792458D5        ; speed of light [km/s]

    root = '92kennicutt'
    
    path = getenv('CATALOGS_DIR')+'/'+root+'/'
    raw = rsex(path+root+'.dat')
    ngalaxy = n_elements(raw)

    ancillary = mrdfits(path+root+'_ancillary_data.fits.gz',1,/silent)
;   ned = mrdfits(path+root+'_ned.fits.gz',1,/silent)
;   photo = mrdfits(path+root+'_ned_photo.fits.gz',1,/silent)

; ---------------------------------------------------------------------------    
; preliminaries    
; ---------------------------------------------------------------------------    
    
;   parse_ned_byname, '92kennicutt.ned', inputnedfile='galaxy_list.txt', $
;     outfile='92kennicutt_ned.fits', nedpath=path, outpath=path
;   parse_ned_photometry, '92kennicutt.photometry.ned', inputnedfile='galaxy_list.txt', $
;     outfile='92kennicutt_ned_photo.fits', nedpath=path, outpath=path

;   basicname = '92kennicutt_ned.fits.gz'
;   photoname = '92kennicutt_ned_photo.fits.gz'
;   distname = '92kennicutt_distances.fits.gz'
;   diamname = '92kennicutt_diameters.fits.gz'
;   
;   write_ancillary_data, ancillary, datapath=outpath, outpath=outpath, $
;     basicname=basicname, photoname=photoname, distname=distname, $
;     diamname=diamname
;   outname = '92kennicutt_ancillary_data.fits'
;   mwrfits, ancillary, outpath+outname, /create
;   spawn, ['gzip -f '+outpath+outname], /sh
    
; ---------------------------------------------------------------------------    

    data = init_cat_linefit(ngalaxy=ngalaxy)
    data = struct_addtags(data,replicate({type: '', aperture: '', scanlen: 0.0, $
      aperwid: 0.0, posangle: 0.0, resolution: '', c41_50: 0.0D},ngalaxy))
    data = struct_addtags(data,struct_trimtags(ancillary,select=[$
      'RC3_B','RC3_B_ERR','RC3_M_B','RC3_M_B_ERR',$
      'RC3_B_LUM','RC3_B_LUM_ERR','RC3_V','RC3_V_ERR','RC3_M_V','RC3_M_V_ERR',$
      'RC3_V_LUM','RC3_V_LUM_ERR','L_FIR_L_B','L_FIR_L_B_ERR']))

; ---------------------------------------------------------------------------    
; properties of interest
; ---------------------------------------------------------------------------    

    data.galaxy     = raw.galaxy
    data.type       = raw.type
    data.aperture   = raw.aperture
    data.scanlen    = raw.scanlen
    data.aperwid    = raw.aperwid
    data.posangle   = raw.posangle
    data.resolution = raw.resolution
    data.c41_50     = raw.c41_50

    data.alt_galaxy = ancillary.ned_galaxy
    data.ra         = ancillary.ra
    data.dec        = ancillary.dec
    data.z_obj      = ancillary.z
    data.distance   = ancillary.distance

;   good = where(photo.b_rc3 gt -900.0)
;   abs = im_absolute_magnitudes('B',photo[good].b_rc3,$
;     data[good].distance,mag_err=photo[good].b_rc3_err)
;   data[good].M_B       = abs.absmag
;   data[good].M_B_err   = abs.absmag_err
;   data[good].B_lum     = abs.lum
;   data[good].B_lum_err = abs.lum_err

; ---------------------------------------------------------------------------    
; EW's; assume 30% uncertainties
; ---------------------------------------------------------------------------    

    good = where(raw.oii_3727_ew gt -900,ngood)
    data[good].oii_3727_ew[0] = raw[good].oii_3727_ew
    data[good].oii_3727_ew[1] = data[good].oii_3727_ew[0]*0.30
    
    good = where(raw.h_beta_ew gt -900,ngood)
    data[good].h_beta_ew[0] = raw[good].h_beta_ew
    data[good].h_beta_ew[1] = data[good].h_beta_ew[0]*0.30
    
    good = where(raw.oiii_5007_ew gt -900,ngood)
    data[good].oiii_5007_ew[0] = raw[good].oiii_5007_ew
    data[good].oiii_5007_ew[1] = data[good].oiii_5007_ew[0]*0.30
    data[good].oiii_4959_ew    = data[good].oiii_5007_ew/3.0

; jm05oct21uofa - adopt a constant [N II]/Ha = 0.5

    good = lindgen(ngalaxy)
    niiha = 0.5
    nii_cor = 1.0 + niiha
    ha_cor = 1.0 + 1.0/niiha
;   good = where((raw.h_alpha_nii_ew gt -900) and (raw.nii_h_alpha gt -900),ngood)
;   nii_cor = 1.0 + raw[good].nii_h_alpha
;   ha_cor = 1.0 + 1.0/raw[good].nii_h_alpha

    data[good].h_alpha_ew[0] = raw[good].h_alpha_nii_ew/nii_cor
    data[good].h_alpha_ew[1] = data[good].h_alpha_ew[0]*0.30
    
    data[good].nii_6584_ew[0] = raw[good].h_alpha_nii_ew/ha_cor
    data[good].nii_6584_ew[1] = data[good].nii_6584_ew[0]*0.30
    data[good].nii_6548_ew    = data[good].nii_6584_ew/3.0
    
; ---------------------------------------------------------------------------
; fluxes
; ---------------------------------------------------------------------------
    
    good = where((raw.oii_3727 gt -900),ngood)
    data[good].oii_3727[0] = raw[good].oii_3727
    data[good].oii_3727[1] = data[good].oii_3727[0]*0.30
    
    good = where((raw.h_beta gt -900),ngood)
    data[good].h_beta[0] = raw[good].h_beta
    data[good].h_beta[1] = data[good].h_beta[0]*0.30
    
    good = where((raw.oiii_5007 gt -900),ngood)
    data[good].oiii_5007[0] = raw[good].oiii_5007
    data[good].oiii_5007[1] = data[good].oiii_5007[0]*0.30
    data[good].oiii_4959    = data[good].oiii_5007/3.0

    good = where((raw.nii_h_alpha gt -900),ngood)
    nii_cor = 1.0 + raw[good].nii_h_alpha
    ha_cor = 1.0 + 1.0/raw[good].nii_h_alpha

    data[good].h_alpha[0] = 1.0/nii_cor
    data[good].h_alpha[1] = data[good].h_alpha[0]*0.30

    data[good].nii_6584[0] = 1.0/ha_cor
    data[good].nii_6584[1] = data[good].nii_6584[0]*0.30
    data[good].nii_6548    = data[good].nii_6584/3.0

;   good = where((raw.oii_3727 gt -900) and (raw.nii_h_alpha gt -900),ngood)
;   nii_cor = 1.0 + raw[good].nii_h_alpha
;   data[good].oii_3727[0] = raw[good].oii_3727*nii_cor
;   data[good].oii_3727[1] = data[good].oii_3727[0]*0.30
;   
;   good = where((raw.h_beta gt -900) and (raw.nii_h_alpha gt -900),ngood)
;   nii_cor = 1.0 + raw[good].nii_h_alpha
;   data[good].h_beta[0] = raw[good].h_beta*nii_cor
;   data[good].h_beta[1] = data[good].h_beta[0]*0.30
;   
;   good = where((raw.oiii_5007 gt -900) and (raw.nii_h_alpha gt -900),ngood)
;   nii_cor = 1.0 + raw[good].nii_h_alpha
;   data[good].oiii_5007[0] = raw[good].oiii_5007*nii_cor
;   data[good].oiii_5007[1] = data[good].oiii_5007[0]*0.30
;   data[good].oiii_4959    = data[good].oiii_5007/3.0
;
;   good = where((raw.nii_h_alpha gt -900),ngood)
;   nii_cor = 1.0 + raw[good].nii_h_alpha
;   data[good].h_alpha[0] = 1.0/nii_cor
;   data[good].h_alpha[1] = data[good].h_alpha[0]*0.30
    
;   nii_ha_ratio = 0.5       ; assume a constant (mean) [N II]/Ha flux ratio
;   nii_cor = 1+nii_ha_ratio
;
;   g = where(c41_50 gt -90.0)    & data[g].c41_50       = c41_50[g]
;   g = where(oii gt -90.0)       & data[g].oii_h_alpha  = alog10(oii[g]*nii_cor)
;   g = where(oiii gt -90.0)      & data[g].oiii_h_alpha = alog10(oiii[g]*nii_cor)
;   g = where(sii gt -90.0)       & data[g].sii_h_alpha  = alog10(sii[g]*nii_cor)
;
;   g = where(oii_ew gt -90.0)    & data[g].oii_ew       = oii_ew[g]
;   g = where(hb_ew gt -90.0)     & data[g].h_beta_ew    = hb_ew[g]
;   g = where(oiii_ew gt -90.0)   & data[g].oiii_ew      = oiii_ew[g]
;   g = where(ha_nii_ew gt -90.0) & data[g].h_alpha_ew   = ha_nii_ew[g]/nii_cor
;   g = where(sii_ew gt -90.0)    & data[g].sii_ew       = sii_ew[g]
    
; ---------------------------------------------------------------------------    
; continuum fluxes
; ---------------------------------------------------------------------------    

    splog, 'Writing '+path+root+'.fits.'
    mwrfits, data, path+root+'.fits', /create
    spawn, ['gzip -f ']+path+root+'.fits', /sh

return
end    
