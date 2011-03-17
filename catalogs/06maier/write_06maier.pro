pro write_06maier, outdata
; jm06apr12uofa - 

    snrcut = 3.0
    Bmsun = k_solar_magnitudes(filterlist='bessell_B.par')
    Bvega2ab = (k_vega2ab(filterlist='bessell_B.par',/kurucz))[0]
    
    root = '06maier'
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
    data.z_obj = raw.z
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

; all (two-sigma) upper limits    
    
    data[good].nii_6584[0] = -raw.nii*1D-17/2.0
    data[good].nii_6584[1] = -3.0
    data[good].nii_6548[0] = data.nii_6584[0] / 3.0 ; assume this is true
    data[good].nii_6548[1] = -3.0

; ---------------------------------------------------------------------------    
; compute the reddening, metallicity, and other miscellaneous
; properties, then write out 
; ---------------------------------------------------------------------------    

    idata = iunred_linedust(data,snrcut=snrcut,/silent,/nopropagate);,cutsig=1.5) ; relax the default CUTSIG
    mz_log12oh, data, idata, data, ohdust, ohnodust, abund, abundnodust, $
      r23branch=r23branch, branchmethod=4L, snrcut=snrcut

    outdata = struct_addtags(data,ohdust)

    lums = compute_linelums(outdata,select_lines=['oii_3727','oiii_5007','h_beta','h_alpha'])
    outdata = struct_addtags(outdata,lums)
    sfrs = compute_sfrs(outdata)

    class = iclassification(outdata,/kauffmann,/silent)
    outdata = struct_addtags(struct_addtags(outdata,sfrs),class)

    splog, 'Writing '+path+root+'.fits.'
    mwrfits, outdata, path+root+'.fits', /create
    spawn, ['gzip -f ']+path+root+'.fits', /sh

; de-reddened data    

    outdata = struct_addtags(idata,ohnodust)

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

    indx = where((outdata.m_b gt -900.0) and (outdata.zstrong_12oh_kk04 gt -900.0),nindx)
    
    print, '##################################################'
    splog, 'Maier et al. (2006): '+string(nindx,format='(I0.0)')+' galaxies.'
    stats = im_stats(outdata[indx].z_obj)
    splog, strtrim(string(stats.min,format='(F12.2)'),2)+'<z<'+strtrim(string(stats.max,format='(F12.2)'),2)+$
      ', <z>='+strtrim(string(stats.mean,format='(F12.2)'),2)
    stats = im_stats(outdata[indx].m_b)
    splog, strtrim(string(stats.min,format='(F12.2)'),2)+'<M_B<'+strtrim(string(stats.max,format='(F12.2)'),2)+$
      ', <M_B>='+strtrim(string(stats.mean,format='(F12.2)'),2)
    print, '##################################################'
    struct_print, struct_trimtags(outdata,select=['GALAXY','M_B','ZSTRONG_12OH_KK04*','R23BRANCH_KK04'])
    print, '##################################################'
       
return
end    
