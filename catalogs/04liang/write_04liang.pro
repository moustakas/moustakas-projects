pro write_04liang, outdata, debug=debug
; jm04nov30uofa

    snrcut = 3.0
    Bmsun = k_solar_magnitudes(filterlist='bessell_B.par')
    Bvega2ab = (k_vega2ab(filterlist='bessell_B.par',/kurucz))[0]
    hahb = 2.86
    
    root = '04liang'
    path = getenv('CATALOGS_DIR')+'/'+root+'/'
    table1 = rsex(path+root+'_table1.dat')
    table2 = rsex(path+root+'_table2.dat')
    liangraw = rsex(path+root+'.dat')

; UDSF13 IS AN AGN    
    
; concatenate tables 1 & 2

    table1_cut = struct_trimtags(table1,except=['CFRS_NAME','I_AB_PHOT','I_AB_SPEC'])
    table2_cut = struct_trimtags(table2,except=['R_AB_PHOT','R_AB_SPEC'])
    table12 = struct_append(table1_cut,table2_cut)

    match, strtrim(liangraw.galaxy,2), strtrim(table12.galaxy,2), indx1, indx2
;   niceprint, liangraw[indx1].galaxy, table12[indx2].galaxy
;   niceprint, liangraw[indx1].galaxy, liangraw[indx1].z, table12[indx2].z

    raw = struct_addtags(liangraw[indx1],struct_trimtags(table12[indx2],except=['GALAXY','Z']))
    ngalaxy = n_elements(raw)

; initialize the output data structure    
    
    data = init_cat_linefit(ngalaxy=ngalaxy)

; ---------------------------------------------------------------------------    
; properties of interest
; ---------------------------------------------------------------------------    

    data.id = 'L04_'+string(lindgen(ngalaxy),format='(I2.2)')
    data.galaxy = raw.galaxy
    data.alt_galaxy = raw.alt_galaxy
    data.ra = raw.ra
    data.dec = raw.dec
    data.z_obj = raw.z
    data.distance = dluminosity(data.z_obj,/Mpc)

    good = where(raw.M_B_AB gt -900,ngood)
    if (ngood ne 0L) then begin
       
       data[good].M_B = raw[good].M_B_AB - Bvega2ab ; convert to Vega
       data[good].M_B_err = 0.1

       lum = 10.0D^(0.4*Bmsun) * 10.0^(-0.4*data[good].M_B) ; [L_sun]
       lum_err = alog(10.0) * 0.4D * lum * data[good].M_B_err ; [L_sun]
       
       data[good].B_lum_err = lum_err/lum/alog(10.0) ; log L_sun
       data[good].B_lum = alog10(lum) ; log L_sun

    endif

    moretags = replicate({lir: -999.0},ngalaxy)
    good = where(raw.lir gt -900.0,ngood)
    if (ngood ne 0L) then moretags[good].lir = raw[good].lir
    data = struct_addtags(data,moretags)

; ---------------------------------------------------------------------------
; fluxes    
; ---------------------------------------------------------------------------

    good = where(raw.oii_3727 gt -900,ngood)
    if (ngood ne 0L) then begin
       data[good].oii_3727 = 1D-17*transpose([ [raw[good].oii_3727], [raw[good].oii_3727_err*raw[good].oii_3727] ])
       data[good].oii_3727[0] = data[good].oii_3727[0]*(1.0+data[good].z_obj)
       data[good].oii_3727[1] = data[good].oii_3727[1]*(1.0+data[good].z_obj)
    endif

    good = where(raw.h_gamma gt -900,ngood)
    if (ngood ne 0L) then begin
       data[good].h_gamma = 1D-17*transpose([ [raw[good].h_gamma], [raw[good].h_gamma_err*raw[good].h_gamma] ])
       data[good].h_gamma[0] = data[good].h_gamma[0]*(1.0+data[good].z_obj)
       data[good].h_gamma[1] = data[good].h_gamma[1]*(1.0+data[good].z_obj)
    endif

    good = where(raw.h_beta gt -900,ngood)
    if (ngood ne 0L) then begin
       data[good].h_beta = 1D-17*transpose([ [raw[good].h_beta], [raw[good].h_beta_err*raw[good].h_beta] ])
       data[good].h_beta[0] = data[good].h_beta[0]*(1.0+data[good].z_obj)
       data[good].h_beta[1] = data[good].h_beta[1]*(1.0+data[good].z_obj)
    endif

    good = where(raw.h_alpha gt -900,ngood)
    if (ngood ne 0L) then begin
       data[good].h_alpha = 1D-17*transpose([ [raw[good].h_alpha], [raw[good].h_alpha_err*raw[good].h_alpha] ])
       data[good].h_alpha[0] = data[good].h_alpha[0]*(1.0+data[good].z_obj)
       data[good].h_alpha[1] = data[good].h_alpha[1]*(1.0+data[good].z_obj)
    endif

; predict H-alpha given H-beta and E(B-V) from Hb/Hg
    
;   ebvgood = where((data.h_beta[1] gt 0.0) and (data.h_gamma[1] gt 0.0),ngood)
;   if (ngood ne 0L) then begin
;      ebv = get_ebv(data[ebvgood].h_beta[0]/data[ebvgood].h_gamma[0],/hbhg)
;      factor = hahb*10^(0.4*ebv*(k_lambda(4861.0,/odonnell)-k_lambda(6563.0,/odonnell)))
;      data[ebvgood].h_alpha[0] = data[ebvgood].h_beta[0]*factor
;      data[ebvgood].h_alpha[1] = data[ebvgood].h_beta[1]*factor
;   endif
    
    good = where(raw.oiii_5007 gt -900,ngood)
    if (ngood ne 0L) then begin
       data[good].oiii_5007 = 1D-17*transpose([ [raw[good].oiii_5007], [raw[good].oiii_5007_err*raw[good].oiii_5007] ])
       data[good].oiii_4959 = data[good].oiii_5007 / 3.0
       data[good].oiii_5007[0] = data[good].oiii_5007[0]*(1.0+data[good].z_obj)
       data[good].oiii_5007[1] = data[good].oiii_5007[1]*(1.0+data[good].z_obj)
       data[good].oiii_4959[0] = data[good].oiii_4959[0]*(1.0+data[good].z_obj)
       data[good].oiii_4959[1] = data[good].oiii_4959[1]*(1.0+data[good].z_obj)
    endif

    good = where(raw.nii_6584 gt -900,ngood)
    if (ngood ne 0L) then begin
       data[good].nii_6584 = 1D-17*transpose([ [raw[good].nii_6584], [raw[good].nii_6584_err*raw[good].nii_6584] ])
       data[good].nii_6548 = data[good].nii_6584 / 3.0
       data[good].nii_6584[0] = data[good].nii_6584[0]*(1.0+data[good].z_obj)
       data[good].nii_6584[1] = data[good].nii_6584[1]*(1.0+data[good].z_obj)
       data[good].nii_6548[0] = data[good].nii_6548[0]*(1.0+data[good].z_obj)
       data[good].nii_6548[1] = data[good].nii_6548[1]*(1.0+data[good].z_obj)
    endif

; ---------------------------------------------------------------------------    
; EW's - assign the same S/N as for the fluxes; upper limits are negative
; ---------------------------------------------------------------------------    

    good = where((data.oii_3727[1] gt 0.0) and (raw.ewoii gt 0.0),ngood)
    if (ngood ne 0L) then begin
;      niceprint, data[good].oii_3727[0]/data[good].oii_3727[1], raw[good].ewoii
       data[good].oii_3727_ew[0] = raw[good].ewoii
       data[good].oii_3727_ew[1] = raw[good].ewoii / (data[good].oii_3727[0]/data[good].oii_3727[1])
    endif
    
    good = where((data.oiii_5007[1] gt 0.0) and (raw.ewoiii gt 0.0),ngood)
    if (ngood ne 0L) then begin
;      niceprint, data[good].oiii_5007[0]/data[good].oiii_5007[1], raw[good].ewoiii
       data[good].oiii_5007_ew[0] = raw[good].ewoiii
       data[good].oiii_5007_ew[1] = raw[good].ewoiii / (data[good].oiii_5007[0]/data[good].oiii_5007[1])
       data[good].oiii_4959_ew = data[good].oiii_5007_ew / 3.0 ; assume this is true
    endif
    
    good = where((data.h_beta[1] gt 0.0) and (raw.ewhb gt 0.0),ngood)
    if (ngood ne 0L) then begin
;      niceprint, data[good].h_beta[0]/data[good].h_beta[1], raw[good].ewhb
       data[good].h_beta_ew[0] = raw[good].ewhb
       data[good].h_beta_ew[1] = raw[good].ewhb / (data[good].h_beta[0]/data[good].h_beta[1])
    endif
    
; ---------------------------------------------------------------------------    
; compute the reddening, metallicity, and other miscellaneous
; properties, then write out 
; ---------------------------------------------------------------------------    

    idata = iunred_linedust(data,snrcut=snrcut,/silent,/nopropagate)
    mz_log12oh, data, idata, data, ohdust, ohnodust, abund, abundnodust, $
      r23branch=r23branch, branchmethod=4L, snrcut=snrcut

;   niceprint, ebv, outdata[ebvgood].ebv_hbhg, outdata[ebvgood].ebv_hahb
    
    outdata = struct_addtags(data,ohdust)

    lums = compute_linelums(outdata,select_lines=['oii_3727','oiii_5007','h_beta','h_alpha'])
    outdata = struct_addtags(outdata,lums)
    sfrs = compute_sfrs(outdata)

    class = iclassification(outdata,/kauffmann)
    outdata = struct_addtags(struct_addtags(outdata,sfrs),class)

    splog, 'Writing '+path+root+'.fits.'
    mwrfits, outdata, path+root+'.fits', /create
    spawn, ['gzip -f ']+path+root+'.fits', /sh

; de-reddened data    

    outdata = struct_addtags(idata,ohnodust)

    lums = compute_linelums(outdata,select_lines=['oii_3727','oiii_5007','h_beta','h_alpha'])
    outdata = struct_addtags(outdata,lums)
    sfrs = compute_sfrs(outdata)

    class = iclassification(outdata,/kauffmann)
    outdata = struct_addtags(struct_addtags(outdata,sfrs),class)

    splog, 'Writing '+path+root+'_nodust.fits.'
    mwrfits, outdata, path+root+'_nodust.fits', /create
    spawn, ['gzip -f ']+path+root+'_nodust.fits', /sh

; ---------------------------------------------------------------------------    
; write out statistics
; ---------------------------------------------------------------------------    

    indx = where((outdata.m_b gt -900.0) and (outdata.zstrong_ew_12oh_kk04 gt -900.0),nindx)
    
    print, '##################################################'
    splog, 'Liang et al. (2004): '+string(nindx,format='(I0.0)')+' galaxies.'
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
