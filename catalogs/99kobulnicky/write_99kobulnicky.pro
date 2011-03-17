pro write_99kobulnicky, outdata, debug=debug
; jm04jun08uofa
; jm04dec31uofa - re-written

    snrcut = 3.0
    Bmsun = k_solar_magnitudes(filterlist='bessell_B.par')
    
    root = '99kobulnicky'
    
    path = getenv('CATALOGS_DIR')+'/'+root+'/'
    raw = rsex(path+root+'.dat')
    ngalaxy = n_elements(raw)
    
    data = init_cat_linefit(ngalaxy=ngalaxy)

; ---------------------------------------------------------------------------    
; properties of interest
; ---------------------------------------------------------------------------    

    data.id = string(lindgen(ngalaxy),format='(I2.2)')
    data.galaxy = raw.galaxy
    data.z_obj = raw.z
    data.distance = dluminosity(data.z_obj,/Mpc)

    data.M_B = raw.M_B - 5.0*alog10(50.0/70.0) ; H_0 = 50 --> 70 km/s/Mpc
    data.M_B_err = 0.1

    lum = 10.0D^(0.4*Bmsun) * 10.0^(-0.4*data.M_B)   ; [L_sun]
    lum_err = alog(10.0) * 0.4D * lum * data.M_B_err ; [L_sun]
       
    data.B_lum_err = lum_err/lum/alog(10.0) ; log L_sun
    data.B_lum = alog10(lum)                ; log L_sun

;   data = struct_addtags(data,replicate({z_12oh_te: -999.0},ngalaxy))
;   data.z_12oh_te = raw.logoh_te

; ---------------------------------------------------------------------------
; manipulate the raw data 
; ---------------------------------------------------------------------------

; undo the reddening correction, using the Seaton extinction law and
; the c(Hb) values in the table

; ---------------------------------------------------------------------------    
; EW's
; ---------------------------------------------------------------------------    

    data.h_beta_ew[0] = raw.h_beta_ew
    data.h_beta_ew[1] = raw.h_beta_ew_err

; ---------------------------------------------------------------------------
; fluxes    
; ---------------------------------------------------------------------------

    raw.f_h_beta = raw.f_h_beta*1D-16
    raw.f_h_beta_err = raw.f_h_beta_err*1D-16

    good = where(raw.h_beta[0] gt -900.0)
    data[good].h_beta[0] = raw[good].h_beta*raw[good].f_h_beta
    data[good].h_beta[1] = raw[good].h_beta_err*raw[good].f_h_beta

    good = where(raw.oii_3727[0] gt -900.0)
    data[good].oii_3727[0] = raw[good].oii_3727*raw[good].f_h_beta
    data[good].oii_3727[1] = raw[good].oii_3727_err*raw[good].f_h_beta

    good = where(raw.oiii_5007[0] gt -900.0)
    data[good].oiii_5007[0] = raw[good].oiii_5007*raw[good].f_h_beta
    data[good].oiii_5007[1] = raw[good].oiii_5007_err*raw[good].f_h_beta
    data[good].oiii_4959 = data[good].oiii_5007 / 3.0D ; assume this is true

; all the fluxes have already been de-reddened

    good = where(raw.h_alpha[0] gt -900.0)
    data[good].h_alpha[0] = raw[good].h_alpha*raw[good].f_h_beta
    data[good].h_alpha[1] = raw[good].h_alpha*raw[good].f_h_beta
    
    good = where(raw.nii_6584[0] gt -900.0)
    data[good].nii_6584[0] = raw[good].nii_6584*raw.f_h_beta
    data[good].nii_6584[1] = raw[good].nii_6584_err*raw.f_h_beta

; ---------------------------------------------------------------------------    
; compute the reddening, metallicity, and other miscellaneous
; properties, then write out; NB: only the de-reddened data are 
; available!
; ---------------------------------------------------------------------------    

    mz_log12oh, data, data, data, ohdust, ohnodust, abund, abundnodust, $
      r23branch=raw.branch, branchmethod=5L, snrcut=snrcut
    outdata = struct_addtags(data,ohdust)

;   adata = im_abundance(data,snrcut=snrcut)
;   outdata = struct_addtags(data,adata)
;   outdata = abundance_catalogs_log12oh(outdata)

    lums = compute_linelums(outdata,select_lines=['oii_3727','oiii_5007','h_beta','h_alpha'])
    outdata = struct_addtags(outdata,lums)
    sfrs = compute_sfrs(outdata)

    class = iclassification(outdata,/kauffmann,/silent)
    outdata = struct_addtags(struct_addtags(outdata,sfrs),class)

    splog, 'Writing '+path+root+'.fits.'
    mwrfits, outdata, path+root+'.fits', /create
    spawn, ['gzip -f ']+path+root+'.fits', /sh

; ---------------------------------------------------------------------------    
; write out statistics
; ---------------------------------------------------------------------------    

    indx = where((outdata.m_b gt -900.0) and (outdata.zstrong_12oh_kk04 gt -900.0),nindx)
    
    print, '##################################################'
    splog, 'Kobulnicky & Zaritsky (1999): '+string(nindx,format='(I0.0)')+' galaxies.'
    stats = im_stats(outdata[indx].z_obj)
    splog, strtrim(string(stats.min,format='(F12.2)'),2)+'<z<'+strtrim(string(stats.max,format='(F12.2)'),2)+$
      ', <z>='+strtrim(string(stats.mean,format='(F12.2)'),2)
    stats = im_stats(outdata[indx].m_b)
    splog, strtrim(string(stats.min,format='(F12.2)'),2)+'<M_B<'+strtrim(string(stats.max,format='(F12.2)'),2)+$
      ', <M_B>='+strtrim(string(stats.mean,format='(F12.2)'),2)
    print, '##################################################'
    struct_print, struct_trimtags(outdata,select=['GALAXY','M_B','ZSTRONG_12OH_KK04*','R23BRANCH_KK04'])
    print, '##################################################'
       
; ---------------------------------------------------------------------------    
; make some plots
; ---------------------------------------------------------------------------    

    if keyword_set(debug) then begin
       
       plot, outdata.m_b, outdata.zstrong_12oh_kk04, xrange=[-23,-15], yrange=[8.0,9.2], $
         ps=4, xsty=3, ysty=3, sym=2
       cc = get_kbrd(1)

    endif
       
return
end    
