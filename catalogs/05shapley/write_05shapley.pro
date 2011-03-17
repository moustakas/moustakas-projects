pro write_05shapley, outdata
; jm04dec02uofa

    snrcut = 3.0
    Bmsun = k_solar_magnitudes(filterlist='bessell_B.par')

    root = '05shapley'
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

    data.M_B = raw.M_B
    data.M_B_err = 0.1

    lum = 10.0D^(0.4*Bmsun) * 10.0^(-0.4*data.M_B)   ; [L_sun]
    lum_err = alog(10.0) * 0.4D * lum * data.M_B_err ; [L_sun]
       
    data.B_lum_err = lum_err/lum/alog(10.0) ; log L_sun
    data.B_lum = alog10(lum)                ; log L_sun

    data.mass     = raw.mass
    data.mass_err = raw.mass_err

; ---------------------------------------------------------------------------
; fluxes    
; ---------------------------------------------------------------------------

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

    good = where(raw.hb gt -900,ngood)
    if (ngood ne 0L) then data[good].h_beta = 1D-17*transpose([ [raw[good].hb], [raw[good].hb_err] ])

    good = where(raw.ha gt -900,ngood)
    if (ngood ne 0L) then data[good].h_alpha = 1D-17*transpose([ [raw[good].ha], [raw[good].ha_err] ])

; ---------------------------------------------------------------------------    
; compute the reddening, metallicity, and other miscellaneous
; properties, then write out 
; ---------------------------------------------------------------------------    

    idata = iunred_linedust(data,snrcut=snrcut,/silent,/nopropagate)
    mz_log12oh, data, idata, data, ohdust, ohnodust, abund, abundnodust, $
      r23branch=r23branch, branchmethod=4L, snrcut=snrcut

    outdata = struct_addtags(struct_addtags(data,ohdust),struct_trimtags(abund,select=['*PETTINI*']))

    lums = compute_linelums(outdata,select_lines=['oii_3727','oiii_5007','h_beta','h_alpha'])
    outdata = struct_addtags(outdata,lums)
    sfrs = compute_sfrs(outdata)

    class = iclassification(outdata,/kauffmann,/silent)
    outdata = struct_addtags(struct_addtags(outdata,sfrs),class)

    splog, 'Writing '+path+root+'.fits.'
    mwrfits, outdata, path+root+'.fits', /create
    spawn, ['gzip -f ']+path+root+'.fits', /sh

; de-reddened data    

    outdata = struct_addtags(struct_addtags(idata,ohnodust),struct_trimtags(abundnodust,select=['*PETTINI*']))

    lums = compute_linelums(outdata,select_lines=['oii_3727','oiii_5007','h_beta','h_alpha'])
    outdata = struct_addtags(outdata,lums)
    sfrs = compute_sfrs(outdata)

    class = iclassification(outdata,/kauffmann,/silent)
    outdata = struct_addtags(struct_addtags(outdata,sfrs),class)

    splog, 'Writing '+path+root+'_nodust.fits.'
    mwrfits, outdata, path+root+'_nodust.fits', /create
    spawn, ['gzip -f ']+path+root+'_nodust.fits', /sh
       
; ---------------------------------------------------------------------------    
; write out statistics
; ---------------------------------------------------------------------------    

    indx = where((outdata.m_b gt -900.0) and (outdata.zstrong_12oh_niiha_pettini gt -900.0),nindx)
    
    print, '##################################################'
    splog, 'Shapley et al. (2005): '+string(nindx,format='(I0.0)')+' galaxies.'
    stats = im_stats(outdata[indx].z_obj)
    splog, strtrim(string(stats.min,format='(F12.2)'),2)+'<z<'+strtrim(string(stats.max,format='(F12.2)'),2)+$
      ', <z>='+strtrim(string(stats.mean,format='(F12.2)'),2)
    stats = im_stats(outdata[indx].m_b)
    splog, strtrim(string(stats.min,format='(F12.2)'),2)+'<M_B<'+strtrim(string(stats.max,format='(F12.2)'),2)+$
      ', <M_B>='+strtrim(string(stats.mean,format='(F12.2)'),2)
    print, '##################################################'
    struct_print, struct_trimtags(outdata,select=['GALAXY','M_B','ZSTRONG_12OH_NIIHA_PETTINI*'])
    print, '##################################################'
       
return
end    
