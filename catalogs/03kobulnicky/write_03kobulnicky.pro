pro write_03kobulnicky, outdata, debug=debug
; jm05jan01uofa - Omega0 = 0.3, OmegaL = 0.7, h100 = 0.7

    snrcut = 3.0
    Bmsun = k_solar_magnitudes(filterlist='bessell_B.par')
    hbeta_abs = 2.0 ; Balmer absorption correction
    hbeta_abs_err = 2.0

    root = '03kobulnicky'
    path = getenv('CATALOGS_DIR')+'/'+root+'/'
    raw = im_read_fmr(path+root+'_table1.dat')

; only keep objects with well-measured emission lines; exclude AGN  

    good = where((strmatch(raw.note,'*B*') eq 0B) and (strmatch(raw.note,'*C*') eq 0B) and $
      (strmatch(raw.note,'*D*') eq 0B),ngalaxy)
    raw = raw[good]

; initialize the output data structure    

    data = init_cat_linefit(ngalaxy=ngalaxy)

; ---------------------------------------------------------------------------    
; store relevant galaxy properties
; ---------------------------------------------------------------------------    

    data.id = 'K03_'+string(raw.id,format='(I2.2)')
    data.galaxy = 'DGSS'+strtrim(raw.gid,2)
    data.ra = string(raw.rah,format='(I2)')+':'+string(raw.ram,format='(I2)')+':'+$
      string(raw.ras,format='(F4.2)')
    data.ra = string(raw.ded,format='(I2)')+':'+string(raw.dem,format='(I2)')+':'+$
      string(raw.des,format='(F4.2)')

    data.z_obj = raw.z
    data.distance = dluminosity(data.z_obj,/Mpc)

    data.M_B = raw.Bmag
    data.M_B_err = 0.05 ; conservative, from Gebhardt et al. (2003)

    lum = 10.0D^(0.4*Bmsun) * 10.0^(-0.4*data.M_B)   ; [L_sun]
    lum_err = alog(10.0) * 0.4D * lum * data.M_B_err ; [L_sun]
       
    data.B_lum_err = lum_err/lum/alog(10.0) ; log L_sun
    data.B_lum = alog10(lum)                ; log L_sun
    
; ---------------------------------------------------------------------------    
; EW's
; ---------------------------------------------------------------------------    

; add a 20% systematic uncertainty (Section 2.1)
    
    good = where(raw.ew3727 gt 0.0)
    data[good].oii_3727_ew[0] = raw[good].ew3727
    data[good].oii_3727_ew[1] = sqrt(raw[good].e_ew3727^2+(0.2*raw[good].ew3727)^2)

; correct for stellar absorption (Section 3.1)    
    
    good = where(raw.ew4861 gt 0.0)
    data[good].h_beta_ew[0] = raw[good].ew4861 + hbeta_abs
    data[good].h_beta_ew[1] = sqrt(raw[good].e_ew4861^2 + hbeta_abs_err^2)

    good = where(raw.ew5007 gt 0.0)
    data[good].oiii_5007_ew[0] = raw[good].ew5007
    data[good].oiii_5007_ew[1] = raw[good].e_ew5007

; for 4959 it doesn't matter which assumption we make    
    
;   good = where(raw.ew4959 gt 0.0)
;   data[good].oiii_4959_ew[0] = raw[good].ew4959
;   data[good].oiii_4959_ew[1] = raw[good].e_ew4959

    data[good].oiii_4959_ew = data[good].oiii_5007_ew / 3.0D
    
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
    splog, 'Kobulnicky et al. (2003): '+string(nindx,format='(I0.0)')+' galaxies.'
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
