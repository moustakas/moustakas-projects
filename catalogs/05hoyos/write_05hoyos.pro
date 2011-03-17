pro write_05hoyos, outdata, debug=debug
; jm06apr12uofa - Omega0 = 0.3, OmegaL = 0.7, h100 = 0.7, Vega mags

    snrcut = 0.0
    Bvega2ab = (k_vega2ab(filterlist='bessell_B.par',/kurucz,/silent))[0]
    B2g01 = +0.0759 ; ^{0.1}g = B+0.0759+0.0620*[(B-V)-0.5870] [AB, Blanton & Roweis]

    root = '05hoyos'
    path = getenv('CATALOGS_DIR')+'/'+root+'/'
    raw = rsex(path+root+'.dat')
    ngalaxy = n_elements(raw)

; initialize the output data structure    

    data = init_cat_linefit(ngalaxy=ngalaxy)

; ---------------------------------------------------------------------------    
; store relevant galaxy properties
; ---------------------------------------------------------------------------    

    data.id = lindgen(ngalaxy)+1L
    data.galaxy = raw.galaxy
    data.ra = 15.0D*im_hms2dec(raw.ra)
    data.dec = im_hms2dec(raw.dec)
    data.z = raw.z

    data.lit_log12oh = raw.log12oh
    data.lit_log12oh_err = raw.log12oh_err
    
    good = where(raw.m_b gt -900,ngood) ; h=0.7, Vega
    if (ngood ne 0L) then begin
       data[good].m_b = raw[good].m_b
       data[good].m_g = data[good].m_b + Bvega2ab + B2g01 ; Vega-->AB; B-->g0.1
    endif
       
; read the DEEP2 and TKRS catalogs and spherematch to get the
; photometry; two are missing from DEEP2 and one from TKRS because of
; incomplete coordinates

    deep = read_deep2(/kcorr)
    spherematch, deep.ra, deep.dec, data.ra, data.dec, 5.0D/3600.0D, m1, m2

    niceprint, raw[m2].galaxy, deep[m1].z, raw[m2].z, $
      deep[m1].ubvri_absmag[1], raw[m2].m_b, deep[m1].mass
    data[m2].mass = deep[m1].mass + im_convert_imf(/from_chabrier)
    data[m2].m_b = deep[m1].ubvri_absmag[1]
    data[m2].m_g = deep[m1].ugriz_absmag[1]
    
    tkrs = read_tkrs(/kcorr)
    spherematch, tkrs.ra, tkrs.dec, data.ra, data.dec, 5.0D/3600.0D, m1, m2
    niceprint, raw[m2].galaxy, tkrs[m1].z, raw[m2].z, $
      tkrs[m1].ubvri_absmag[1], raw[m2].m_b, tkrs[m1].mass
    data[m2].mass = tkrs[m1].mass + im_convert_imf(/from_chabrier)
    data[m2].m_b = tkrs[m1].ubvri_absmag[1]
    data[m2].m_g = tkrs[m1].ugriz_absmag[1]

; ---------------------------------------------------------------------------    
; EW's
; ---------------------------------------------------------------------------    

    good = where(raw.ewoii gt 0.0)
    data[good].oii_3727_ew[0] = raw[good].ewoii
    data[good].oii_3727_ew[1] = raw[good].ewoii_err

    good = where(raw.ewhb gt 0.0)
    data[good].h_beta_ew[0] = raw[good].ewhb
    data[good].h_beta_ew[1] = raw[good].ewhb_err

    good = where(raw.ewoiii gt 0.0)
    data[good].oiii_5007_ew[0] = raw[good].ewoiii
    data[good].oiii_5007_ew[1] = raw[good].ewoiii_err

    data[good].oiii_4959_ew = data[good].oiii_5007_ew / 3.0 ; assume this is true
    
; ---------------------------------------------------------------------------    
; assign R23 branches
; ---------------------------------------------------------------------------    

    splog, 'Assigning R23 branches'
    data.r23branch = raw.branch
    
; ---------------------------------------------------------------------------    
; compute the reddening, metallicity, and other miscellaneous
; properties, then write out 
; ---------------------------------------------------------------------------    

    mz_log12oh, data, nmonte=500L, snrcut=0.0, final_ohdust=log12oh, $
      logfile=path+'mz_log12oh.05hoyos.log'
    outdata = struct_addtags(data,log12oh)
    
    splog, 'Writing '+path+root+'.fits'
    mwrfits, outdata, path+root+'.fits', /create

; ---------------------------------------------------------------------------    
; write out statistics
; ---------------------------------------------------------------------------    

    indx = where((outdata.m_b gt -900.0) and (outdata.zstrong_ew_alpha_unity_12oh_kk04 gt -900.0),nindx)
    
    print, '##################################################'
    splog, 'Hoyos et al. (2005): '+string(nindx,format='(I0.0)')+' galaxies.'
    stats = im_stats(outdata[indx].z)
    splog, strtrim(string(stats.min,format='(F12.2)'),2)+'<z<'+strtrim(string(stats.max,format='(F12.2)'),2)+$
      ', <z>='+strtrim(string(stats.mean,format='(F12.2)'),2)
    stats = im_stats(outdata[indx].m_b)
    splog, strtrim(string(stats.min,format='(F12.2)'),2)+'<M_B<'+strtrim(string(stats.max,format='(F12.2)'),2)+$
      ', <M_B>='+strtrim(string(stats.mean,format='(F12.2)'),2)
    print, '##################################################'
    struct_print, struct_trimtags(outdata,select=['GALAXY','M_B',$
      'ZSTRONG_EW_ALPHA_UNITY_12OH_KK04*','R23BRANCH_EW_ALPHA_UNITY_KK04'])
    print, '##################################################'
       
; ---------------------------------------------------------------------------    
; make some plots
; ---------------------------------------------------------------------------    

    if keyword_set(debug) then begin
       
       djs_plot, outdata.m_g, outdata.zstrong_ew_alpha_unity_12oh_kk04, xrange=[-23.5,-17], $
         yrange=[7.0,9.5], ps=4, xsty=3, ysty=3, sym=2, charsize=2.0, charthick=2.0, $
         xtitle='M_{0.1g}', ytitle='12 + log (O/H)'
       djs_oplot, outdata.m_g, outdata.lit_log12oh, ps=7, sym=2.0, color='red'
       cc = get_kbrd(1)

    endif
       
return
end    
