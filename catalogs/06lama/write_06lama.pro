pro write_06lama, outdata, debug=debug
; jm06apr11uofa - 

    snrcut = 3.0
    Bmsun = k_solar_magnitudes(filterlist='bessell_B.par')
    Bvega2ab = (k_vega2ab(filterlist='bessell_B.par',/kurucz))[0]
    
    root = '06lama'
    path = getenv('CATALOGS_DIR')+'/'+root+'/'

    objlist = mrdfits(path+root+'_objlist.fits',1,/silent)
    class = mrdfits(path+root+'_class.fits',1,/silent)
    phot = mrdfits(path+root+'_phot.fits',1,/silent)
    forb = mrdfits(path+root+'_forbidden.fits',1,/silent)
    balmer = mrdfits(path+root+'_balmer.fits',1,/silent)
    lamaraw = rsex(path+root+'.dat')

; only retain the 118 galaxies that appear in the second paper of this
; series (basically, just the star-forming and candidate star-forming
; galaxies)

    match, strtrim(lamaraw.galaxy,2), 'LCL05_'+string(objlist.lcl05,format='(I3.3)'), indx1, indx2
;   niceprint, lamaraw[indx1].m_b, phot[indx2].imag+phot[indx2]._v_i_0+phot[indx2]._b_v_0

    raw = struct_addtags(struct_addtags(lamaraw[indx1],struct_trimtags(forb[indx2],except=['RECNO','LCL05'])),$
      struct_trimtags(balmer[indx2],except=['RECNO','LCL05']))
    ngalaxy = n_elements(raw)

; initialize the output data structure    

    data = init_cat_linefit(ngalaxy=ngalaxy)
    
; ---------------------------------------------------------------------------    
; store relevant galaxy properties
; ---------------------------------------------------------------------------    

    data.id = string(lindgen(ngalaxy),format='(I3.3)')
    data.galaxy = strtrim(raw.galaxy,2)
    data.alt_galaxy = strtrim(raw.alt_galaxy,2)
    data.z_obj = raw.z
    data.distance = dluminosity(data.z_obj,/Mpc)

    good = where(raw.M_B_AB gt -900,ngood)
    if (ngood ne 0L) then begin
       
       data[good].M_B = raw[good].M_B_AB - 5.0*alog10(71.0/70.0) - Bvega2ab ; convert to Vega
       data[good].M_B_err = 0.1

       lum = 10.0D^(0.4*Bmsun) * 10.0^(-0.4*data[good].M_B) ; [L_sun]
       lum_err = alog(10.0) * 0.4D * lum * data[good].M_B_err ; [L_sun]
       
       data[good].B_lum_err = lum_err/lum/alog(10.0) ; log L_sun
       data[good].B_lum = alog10(lum) ; log L_sun

    endif
       
; ---------------------------------------------------------------------------
; fluxes    
; ---------------------------------------------------------------------------

    good = where(raw.fhb gt 0.0)
    data[good].h_beta[0] = raw[good].fhb*1D-17
    data[good].h_beta[1] = raw[good].e_fhb*1D-17

    good = where((raw.fha gt 0.0) and (data.h_beta[1] gt 0.0))
    data[good].h_alpha[0] = raw[good].fha*data[good].h_beta[0]
    data[good].h_alpha[1] = raw[good].e_fha*data[good].h_beta[1]
    
    good = where((raw.fhg gt 0.0) and (data.h_beta[1] gt 0.0))
    data[good].h_gamma[0] = raw[good].fhg*data[good].h_beta[0]
    data[good].h_gamma[1] = raw[good].e_fhg*data[good].h_beta[1]
    
    good = where((raw.fo2 gt 0.0) and (data.h_beta[1] gt 0.0))
    data[good].oii_3727[0] = raw[good].fo2*data[good].h_beta[0]
    data[good].oii_3727[1] = raw[good].e_fo2*data[good].h_beta[1]

    good = where((raw.fo3b gt 0.0) and (data.h_beta[1] gt 0.0))
    data[good].oiii_5007[0] = raw[good].fo3b*data[good].h_beta[0]
    data[good].oiii_5007[1] = raw[good].e_fo3b*data[good].h_beta[1]
    data[good].oiii_4959 = data[good].oiii_5007 / 3.0 ; assume this is true

    good = where((raw.fn2 gt 0.0) and (data.h_beta[1] gt 0.0))
    data[good].nii_6584[0] = raw[good].fn2*data[good].h_beta[0]
    data[good].nii_6584[1] = raw[good].e_fn2*data[good].h_beta[1]
    data[good].nii_6548 = data[good].nii_6584 / 3.0 ; assume this is true

; ---------------------------------------------------------------------------    
; EW's
; ---------------------------------------------------------------------------    

    good = where(raw.ewo2 ne 0.0 and raw.e_ewo2 ne 0.0)
    data[good].oii_3727_ew[0] = -raw[good].ewo2
    data[good].oii_3727_ew[1] = raw[good].e_ewo2

    good = where(raw.ewo3b ne 0.0 and raw.e_ewo3b ne 0.0)
    data[good].oiii_5007_ew[0] = -raw[good].ewo3b
    data[good].oiii_5007_ew[1] = raw[good].e_ewo3b
    data[good].oiii_4959_ew = data[good].oiii_5007_ew / 3.0 ; assume this is true

    good = where(raw.ewn2 ne 0.0 and raw.e_ewn2 ne 0.0)
    data[good].nii_6584_ew[0] = -raw[good].ewn2
    data[good].nii_6584_ew[1] = raw[good].e_ewn2
    data[good].nii_6548_ew = data[good].nii_6584_ew / 3.0 ; assume this is true

    good = where(raw.ewha ne 0.0 and raw.e_ewha ne 0.0)
    data[good].h_alpha_ew[0] = -raw[good].ewha
    data[good].h_alpha_ew[1] = raw[good].e_ewha

    good = where(raw.ewhb ne 0.0 and raw.e_ewhb ne 0.0)
    data[good].h_beta_ew[0] = -raw[good].ewhb
    data[good].h_beta_ew[1] = raw[good].e_ewhb

    good = where(raw.ewhg ne 0.0 and raw.e_ewhg ne 0.0)
    data[good].h_gamma_ew[0] = -raw[good].ewhg
    data[good].h_gamma_ew[1] = raw[good].e_ewhg

; ---------------------------------------------------------------------------    
; compute the reddening, metallicity, and other miscellaneous
; properties, then write out 
; ---------------------------------------------------------------------------    

    idata = iunred_linedust(data,snrcut=snrcut,/silent,/nopropagate)
    mz_log12oh, data, idata, data, ohdust, ohnodust, abund, abundnodust, $
      r23branch=r23branch, branchmethod=4L, snrcut=snrcut

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
    splog, 'Lamareille et al. (2006): '+string(nindx,format='(I0.0)')+' galaxies.'
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
