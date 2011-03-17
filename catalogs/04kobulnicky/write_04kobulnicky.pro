pro write_04kobulnicky, outdata, debug=debug
; jm05jan01uofa - Omega0 = 0.3, OmegaL = 0.7, h100 = 0.7

; taking S/N=3 eliminates #25 and #164, but just barely, so could
; possibly lower the S/N cut    
    
    snrcut = 0.0 ; 2.9

    Bmsun = k_solar_magnitudes(filterlist='bessell_B.par')
    hbeta_abs = 2.0 ; Balmer absorption correction
    hbeta_abs_err = 2.0

    oratio = 2.984
    ocor = 1+(1.0/oratio)
    niihacut = -1.0 ; -1.05 ; NOTE!

    root = '04kobulnicky'
    path = getenv('CATALOGS_DIR')+'/'+root+'/'
    raw = im_read_fmr(path+root+'_table1.dat')

; meaning of the codes:
;    0 = oxygen abundance computed using EWR_23_ = (W_[O II]_ (3727) + 
;        1.3*W_[O III]_ (5007)/W_H{beta}_;
;    1 = object is probably on the turnaround region of the R_23_ - O/H 
;        calibration or possibly on the lower branch based on the 
;        [N II] (6584)/H{alpha} ratio, so lower branch formula of McGaugh (1991)
;        is applied and averaged with the upper branch estimate of KD03 to 
;        arrive at the oxygen abundance in O/H-avg;
;    2 = oxygen abundance computed using EWR_23_ = (W_[O II]_ (3727) + 
;        4 * W_[O III]_ (4959)/W_H{beta}_;
;    3 = W_[O III]_ (5007) is measured with signal-to-noise < 3:1 so a 3{sigma} 
;        lower limit on the oxygen abundance is computed using 
;        EWR_23_ = (W_[O II]_ (3727) + 1.3 * W_[O III]_ (5007))/W_H{beta}_;
    
; remove low S/N (code=3) objects; there are only two objects with
; code=1 (both on the lower branch), and the code=2 objects should be
; fine

    ngalaxy = n_elements(raw)
;   good = where(raw.code eq 0L,ngalaxy)
    good = where(raw.code ne 3L,ngalaxy)
    raw = raw[good]

; initialize the output data structure    

    data = init_cat_linefit(ngalaxy=ngalaxy)
    moretags = replicate({m_b_kk04: -999.0},ngalaxy)
    data = struct_addtags(data,moretags)
    
; ---------------------------------------------------------------------------    
; store relevant galaxy properties
; ---------------------------------------------------------------------------    

    data.id = raw.num
    data.galaxy = raw.goods
    data.z = raw.z

;   allra = strmid(strtrim(raw.goods,2),1,9)
;   alldec = strmid(strtrim(raw.goods,2),10,9)
;   data.ra = 15.0D*im_hms2dec(strmid(allra,0,2)+':'+$
;     strmid(allra,2,2)+':'+strmid(allra,4,5))
;   data.dec = im_hms2dec(strmid(allra,0,3)+':'+$
;     strmid(allra,3,2)+':'+strmid(allra,5,5))
    
    data.lit_log12oh = raw.o_h
    data.lit_log12oh_err = raw.e_o_h
    
    data.m_b_kk04 = raw.Bmag ; Vega, h=0.7

; read the TKRS catalog and spherematch to get the photometry

    tkrs = read_tkrs(/kcorr)

;   spherematch, tkrs.ra, tkrs.dec, data.ra, data.dec, 10.0D/3600.0D, m1, m2
    match, tkrs.id, raw.id, m1, m2
    niceprint, raw[m2].goods, tkrs[m1].goods, tkrs[m1].z, raw[m2].z, $
      tkrs[m1].ubvri_absmag[1], raw[m2].bmag, tkrs[m1].ubvri_absmag[0]-tkrs[m1].ubvri_absmag[1], $
      raw[m2].u_b, tkrs[m1].mass, raw[m2].imag, tkrs[m1].imag_magauto

    data[m2].ra = tkrs[m1].ra
    data[m2].dec = tkrs[m1].dec
    data[m2].mass = tkrs[m1].mass + im_convert_imf(/from_chabrier)
    data[m2].m_b = tkrs[m1].ubvri_absmag[1]
    data[m2].m_g = tkrs[m1].ugriz_absmag[1]
    data[m2].m_r = tkrs[m1].ugriz_absmag[2]
    
; ---------------------------------------------------------------------------    
; EW's
; ---------------------------------------------------------------------------    

    bad = where(raw.e_ewoii eq 0.0,nbad) ; one of these measurements has an uncertainty of zero!
    if (nbad ne 0L) then raw[bad].e_ewoii = 0.01*raw[bad].ewoii ; 1% uncertainty

    bad = where(raw.e_ewoiii eq 0.0,nbad) ; and lots of these!
    if (nbad ne 0L) then raw[bad].e_ewoiii = 0.01*raw[bad].ewoiii ; 1% uncertainty
    
    good = where(raw.ewoii gt 0.0)
    data[good].oii_3727_ew[0] = raw[good].ewoii
    data[good].oii_3727_ew[1] = raw[good].e_ewoii

    good = where(raw.ewhb gt 0.0)
    data[good].h_beta_ew[0] = raw[good].ewhb + hbeta_abs ; correct for stellar absorption
    data[good].h_beta_ew[1] = sqrt(raw[good].e_ewhb^2 + hbeta_abs_err^2)

    good = where(raw.ewoiii gt 0.0)
    data[good].oiii_5007_ew[0] = raw[good].ewoiii / ocor   ; remove the effect of 4959 from their table
    data[good].oiii_5007_ew[1] = raw[good].e_ewoiii / ocor

    data[good].oiii_4959_ew = data[good].oiii_5007_ew / oratio ; assume this is true
    
; ---------------------------------------------------------------------------    
; assign R23 branches
; ---------------------------------------------------------------------------    

    splog, 'Assigning R23 branches'
    data.r23branch = 'U'
    lo = where(raw.code eq 1,nlo)
    if (nlo ne 0L) then data[lo].r23branch = 'L'
    
; ---------------------------------------------------------------------------    
; compute abundances, assigned branches, other miscellaneous
; properties, and write out 
; ---------------------------------------------------------------------------    
    
    ages_mz_log12oh, data, nmonte=500L, snrcut=0.0, final_ohdust=log12oh, $
      logfile=path+'mz_log12oh.04kobulnicky.log'
    outdata = struct_addtags(data,log12oh)

; the biggest outlier is 3272, which for some reason KK04 adopts the
; *average* of the M91 and KD02 abundances; the M91 abundance agrees
; with my number, ~8.67    
    
    splog, 'Writing '+path+root+'.fits'
    mwrfits, outdata, path+root+'.fits', /create

; I don't want to store either [NII] or H-alpha in the final
; structure, since only the ratio is given in the paper; so place a
; handful of objects on the lower branch "by hand", here; use upper
; limits, since the sign goes in the right direction; also, some
; objects are in the turn-around region ('T'), and placing them on the
; lower branch is wrong (e.g., #41), so only re-assign objects to the
; upper branch

; update (jm06jun18uofa) - just do a straight NII/Ha cut   
    
; update (jm06nov21nyu) - updated to use MZ_LOG12OH; place 3792 on the
;                         lower branch based on his luminosity; place
;                         2336 on the upper branch because he is an
;                         outlier on the lower branch
;   
;   ewbranch = strarr(ngalaxy)
;   ewbranch[*] = 'U'
;
;   good = where((raw.ewha_nii gt -900.0),ngood)
;   if (ngood ne 0L) then begin
;      lo = where((-alog10(raw[good].ewha_nii) lt niihacut),nlo)
;      if (nlo ne 0L) then ewbranch[good[lo]] = 'L'
;   endif
;   match = where(raw.id eq 3792,nmatch)
;   if (nmatch ne 0L) then ewbranch[match] = 'L'
;   match = where(raw.id eq 2336,nmatch)
;   if (nmatch ne 0L) then ewbranch[match] = 'U'
;   niceprint, raw[good].num, -alog10(raw[good].ewha_nii), raw[good].o_h_kd
;
;   mz_log12oh, data, data, data, ohdust, ohnodust, abund, abundnodust, $
;     r23branch=ewbranch, branchmethod=5L, snrcut=snrcut
;
;   outdata = struct_addtags(data,ohdust)
    
;   adata = im_abundance(data,snrcut=snrcut)
;   outdata = struct_addtags(data,adata)
;   outdata = abundance_catalogs_log12oh(outdata,ewbranch=ewbranch,branch=ewbranch)

; OLD code
    
;   u = where(outdata.r23branch_kk04_ew eq 'U',nu)
;   l = where(outdata.r23branch_kk04_ew eq 'L',nl)
;   t = where(outdata.r23branch_kk04_ew eq 'T',nt)
; ---
;   ploterror, outdata.zstrong_ew_r23, outdata.zstrong_ew_12oh_kk04_r23_upper, $
;     outdata.zstrong_ew_r23_err, outdata.zstrong_ew_12oh_kk04_r23_upper_err, ps=3, xsty=3, $
;     ysty=3, xr=[0.6,1.0], yr=[8.2,8.7] ; xr=[0,1.1], yr=[7.4,9.3]
;   oploterror, outdata.zstrong_ew_r23, outdata.zstrong_ew_12oh_kk04_r23_lower, $
;     outdata.zstrong_ew_r23_err, outdata.zstrong_ew_12oh_kk04_r23_lower_err, ps=3
; ---
;   plot, outdata.zstrong_ew_r23, outdata.zstrong_ew_12oh_kk04_r23_upper, ps=4, xsty=3, $
;     ysty=3, xr=[0.6,1.0], yr=[8.2,8.7];, xr=[0,1.1], yr=[7.4,9.3]
;   oplot, outdata.zstrong_ew_r23, outdata.zstrong_ew_12oh_kk04_r23_lower, ps=4
; ---
;   plot, outdata[u].zstrong_ew_r23, outdata[u].zstrong_ew_12oh_kk04, ps=4, xsty=3, ysty=3, $
;     xr=[0,1.1], yr=[7.4,9.3]
;   oplot, outdata[u].zstrong_ew_r23, outdata[u].zstrong_ew_12oh_kk04_r23_upper, ps=3
;   oplot, outdata[u].zstrong_ew_r23, outdata[u].zstrong_ew_12oh_kk04_r23_lower, ps=3
; ---
;   if (nl ne 0L) then oploterror, outdata[l].zstrong_ew_r23, outdata[l].zstrong_ew_12oh_kk04, $
;     outdata[l].zstrong_ew_12oh_kk04_err, ps=7, color=djs_icolor('red'), errcolor=djs_icolor('red')
;   if (nt ne 0L) then oploterror, outdata[t].zstrong_ew_r23, outdata[t].zstrong_ew_12oh_kk04, $
;     outdata[t].zstrong_ew_12oh_kk04_err, ps=5, color=djs_icolor('dark green'), errcolor=djs_icolor('dark green')
;   niceprint, outdata.zstrong_ew_r23, outdata.zstrong_ew_12oh_kk04, outdata.r23branch_kk04_ew
    
; ---------------------------------------------------------------------------    
; write out statistics
; ---------------------------------------------------------------------------    

    indx = where((outdata.m_b gt -900.0) and (outdata.zstrong_ew_alpha_unity_12oh_kk04 gt -900.0),nindx)
    
    print, '##################################################'
    splog, 'Kobulnicky & Kewley (2004): '+string(nindx,format='(I0.0)')+' galaxies.'
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
       
       plot, outdata.zstrong_ew_alpha_unity_12oh_kk04, outdata.lit_log12oh, $
         xr=[8.2,9.5], yr=[8.2,9.5], ps=4, charsize=2.0, charthick=2.0, sym=2
       cc = get_kbrd(1)

       plot, outdata.m_b, outdata.zstrong_ew_alpha_unity_12oh_kk04, $
         xrange=[-23,-15], yrange=[8.0,9.2], $
         ps=4, xsty=3, ysty=3, sym=2, charsize=2.0, charthick=2.0
       cc = get_kbrd(1)

    endif
       
return
end    
