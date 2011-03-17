pro write_06mouhcine, outdata, debug=debug
; jm06apr12uofa - Omega0 = 0.3, OmegaL = 0.7, h100 = 0.7, Vega mags

    snrcut = 3.0

    Bmsun = k_solar_magnitudes(filterlist='bessell_B.par')

    root = '06mouhcine'
    path = getenv('CATALOGS_DIR')+'/'+root+'/'
    raw = rsex(path+root+'.dat')
    ngalaxy = n_elements(raw)

; initialize the output data structure    

    data = init_cat_linefit(ngalaxy=ngalaxy)
    
; ---------------------------------------------------------------------------    
; store relevant galaxy properties
; ---------------------------------------------------------------------------    

    data.id = 'MBAN06_'+string(raw.id,format='(I3.3)')
    data.galaxy = data.id
    data.z_obj = raw.z
    data.distance = dluminosity(data.z_obj,/Mpc)

    good = where(raw.M_B gt -900,ngood)
    if (ngood ne 0L) then begin
       
       data[good].M_B = raw[good].M_B
       data[good].M_B_err = raw[good].M_B_err

       lum = 10.0D^(0.4*Bmsun) * 10.0^(-0.4*data[good].M_B) ; [L_sun]
       lum_err = alog(10.0) * 0.4D * lum * data[good].M_B_err ; [L_sun]
       
       data[good].B_lum_err = lum_err/lum/alog(10.0) ; log L_sun
       data[good].B_lum = alog10(lum) ; log L_sun

    endif
       
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
; compute the reddening, metallicity, and other miscellaneous
; properties, then write out 
; ---------------------------------------------------------------------------    

    mz_log12oh, data, data, data, ohdust, ohnodust, abund, abundnodust, $
      r23branch=r23branch, branchmethod=4L, snrcut=snrcut

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

    indx = where((outdata.m_b gt -900.0) and (outdata.zstrong_ew_12oh_kk04 gt -900.0),nindx)
    
    print, '##################################################'
    splog, 'Mouhcine et al. (2006): '+string(nindx,format='(I0.0)')+' galaxies.'
    stats = im_stats(outdata[indx].z_obj)
    splog, strtrim(string(stats.min,format='(F12.2)'),2)+'<z<'+strtrim(string(stats.max,format='(F12.2)'),2)+$
      ', <z>='+strtrim(string(stats.mean,format='(F12.2)'),2)
    stats = im_stats(outdata[indx].m_b)
    splog, strtrim(string(stats.min,format='(F12.2)'),2)+'<M_B<'+strtrim(string(stats.max,format='(F12.2)'),2)+$
      ', <M_B>='+strtrim(string(stats.mean,format='(F12.2)'),2)
    print, '##################################################'
    struct_print, struct_trimtags(outdata,select=['GALAXY','M_B','ZSTRONG_EW_12OH_KK04*','R23BRANCH_EW_KK04'])
    print, '##################################################'
       
; ---------------------------------------------------------------------------    
; make some plots
; ---------------------------------------------------------------------------    

    if keyword_set(debug) then begin
       
       plot, outdata.m_b, outdata.zstrong_ew_12oh_kk04, xrange=[-23,-15], yrange=[8.0,9.2], $
         ps=4, xsty=3, ysty=3, sym=2
       cc = get_kbrd(1)

    endif
       
return
end    
