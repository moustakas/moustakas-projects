pro write_03lilly, outdata, debug=debug, merge_with_maier=merge_with_maier
; jm05jan01uofa - the optical data from this paper has been merged
;                 with the near-infrared observations by Maier et
;                 al. (2005); note that Lilly et al. did not use the
;                 cosmology advertised; hence, their magnitudes are
;                 too bright by 0.45 mag
; jm08apr23nyu - actually, treat the Maier et al. (2005) observations
;                as separate/independent, since the line-fluxes change 

    snrcut = 0.0
    Bmsun = k_solar_magnitudes(filterlist='bessell_B.par')
    Bvega2ab = (k_vega2ab(filterlist='bessell_B.par',/kurucz))[0]
    
    root = '03lilly'
    path = getenv('CATALOGS_DIR')+'/'+root+'/'
    lillyraw = rsex(path+root+'.dat')
;   lillyraw = lillyraw[where(lillyraw.hb_ew gt 5.0)] ; apply EW(Hb) cut as in Lilly et al. 
    
    ngalaxy = n_elements(lillyraw)

    data = init_cat_linefit(ngalaxy=ngalaxy)
    moretags = replicate({m_b_lilly: -999.0, mass_glaze: -999.0, $
      mass_glaze_err: -999.0},ngalaxy)
    data = struct_addtags(data,moretags)

; K-corrections and stellar masses; only trust the results if there is
; photometry in three or more bandpasses; convert to Salpeter 

;   mass = mrdfits(path+'mass_03lilly_kcorr.fits.gz',1,/silent)

    kcorr = mrdfits(path+'03lilly_kcorr.fits.gz',1,/silent)
    match, kcorr.galaxy, lillyraw.galaxy, indx1, indx2
    data = struct_addtags(struct_trimtags(data,except='MASS'),$
      struct_trimtags(kcorr[indx1],except=['GALAXY','Z']))

    nband = total((data.abmaggies gt 0.0),1)
    good = where((data.mass gt -900.0) and (nband ge 3.0))
    data[good].mass = data[good].mass + im_convert_imf(/from_chabrier)
    data[good].m_g = data[good].ugriz_absmag[1]
    data[good].m_r = data[good].ugriz_absmag[2]
    data[good].m_b = data[good].ubvri_absmag[1]

; --------------------------------------------------
; code to match the lilly and maier catalogs
; --------------------------------------------------

    merge_with_maier = 0

    if keyword_set(merge_with_maier) then begin
       
       maierraw = rsex(getenv('CATALOGS_DIR')+'/05maier/05maier.dat')
       bigmaier = im_empty_structure(maierraw,ncopies=n_elements(lillyraw))
       match, strtrim(lillyraw.galaxy,2), strtrim(maierraw.galaxy,2), indx1, indx2
       
       bigmaier[indx1] = maierraw[indx2]

; merge the two tables
       
;      niceprint, bigmaier[indx1].galaxy, lillyraw[indx1].galaxy, bigmaier[indx1].redshift, lillyraw[indx1].z
;      niceprint, bigmaier[indx1].galaxy, lillyraw[indx1].m_b_ab, bigmaier[indx1].m_b_ab

       splog, 'H-beta:'
       niceprint, bigmaier[indx1].galaxy, lillyraw[indx1].hb, bigmaier[indx1].hb, $
         lillyraw[indx1].hb/bigmaier[indx1].hb
       print
       
       splog, '[O II]:'
       niceprint, bigmaier[indx1].galaxy, lillyraw[indx1].oii, bigmaier[indx1].oii, $
         lillyraw[indx1].oii/bigmaier[indx1].oii
       print
       
       splog, '[O III]:'
       niceprint, bigmaier[indx1].galaxy, lillyraw[indx1].oiii, bigmaier[indx1].oiii, $
         lillyraw[indx1].oiii/bigmaier[indx1].oiii
       print

       maierselect = ['HA','HA_ERR','NII','NII_ERR','HA_EW']
       raw = struct_addtags(lillyraw,struct_trimtags(bigmaier,select=maierselect))
;      raw = struct_addtags(lillyraw,im_struct_trimtags(bigmaier,select=tag_names(bigmaier),newtags='MAIER_'+tag_names(bigmaier)))

; the Maier fluxes are presumably superior; not sure whether to change
; the Lilly EWs 

;      raw[indx1].oii      = bigmaier[indx1].oii
;      raw[indx1].oii_err  = bigmaier[indx1].oiii_err
;      raw[indx1].hb       = bigmaier[indx1].hb
;      raw[indx1].hb_err   = bigmaier[indx1].hb_err
;      raw[indx1].oiii     = bigmaier[indx1].oiii
;      raw[indx1].oiii_err = bigmaier[indx1].oiii_err
    
    endif else raw = lillyraw

; ---------------------------------------------------------------------------    
; store relevant galaxy properties
; ---------------------------------------------------------------------------    

    data.id = lindgen(ngalaxy)+1L
    data.galaxy = raw.galaxy
    data.z = raw.z

    good = where(raw.mass gt -900.0)
    data[good].mass_glaze = raw[good].mass + im_convert_imf(/from_baldry03) ; Salpeter
    data[good].mass_glaze_err = raw[good].mass_err

    good = where(raw.m_b_ab gt -900.0)
    data[good].m_b_lilly = raw[good].m_b_ab - Bvega2ab + 0.45 ; AB-->Vega + cosmology shift!

; ---------------------------------------------------------------------------
; fluxes    
; ---------------------------------------------------------------------------
    
    good = where(raw.oii gt 0.0)
    data[good].oii_3727[0] = raw[good].oii*1D-17
    data[good].oii_3727[1] = raw[good].oii_err*1D-17

    good = where(raw.hb gt 0.0)
    data[good].h_beta[0] = raw[good].hb*1D-17 ; already corrected for stellar absorption
    data[good].h_beta[1] = raw[good].hb_err*1D-17

;   data[good].h_alpha = HaHb*data[good].h_beta ; E(B-V)=0
;   data[good].h_alpha[0] = 3.355*raw[good].hb*1D-17 ; this gives E(B-V) = 0.15
;   data[good].h_alpha[1] = 3.355*raw[good].hb_err*1D-17

    good = where(raw.oiii gt 0.0)
    data[good].oiii_5007[0] = raw[good].oiii*1D-17
    data[good].oiii_5007[1] = raw[good].oiii_err*1D-17

    data[good].oiii_4959 = data[good].oiii_5007 / 3.0 ; assume this is true

    if keyword_set(merge_with_maier) then begin
       good = where(raw.ha gt 0.0)
       data[good].h_alpha[0] = raw[good].ha*1D-17
       data[good].h_alpha[1] = raw[good].ha_err*1D-17
       
       good = where(raw.nii gt 0.0)
       data[good].nii_6584[0] = raw[good].nii*1D-17
       data[good].nii_6584[1] = raw[good].nii_err*1D-17

       data[good].nii_6548 = data[good].nii_6584 / 3.0 ; assume this is true
    endif

; ---------------------------------------------------------------------------    
; continuum fluxes; assume the continuum error is the same as the
; flux; also assume that H-beta and [O III] have the same continuum,
; and that H-alpha and [N II] have the same continuum
; ---------------------------------------------------------------------------    

    good = where(raw.oii gt 0.0 and raw.oii_ew gt 0.0)
    data[good].oii_3727_continuum[0] = raw[good].oii/raw[good].oii_ew*1D-17
    data[good].oii_3727_continuum[1] = raw[good].oii_err*1D-17

    good = where(raw.hb gt 0.0 and raw.hb_ew gt 0.0)
    data[good].h_beta_continuum[0] = raw[good].hb/raw[good].hb_ew*1D-17
    data[good].h_beta_continuum[1] = raw[good].hb_err*1D-17
    data[good].oiii_4959_continuum = data[good].h_beta_continuum
    data[good].oiii_5007_continuum = data[good].h_beta_continuum

    if keyword_set(merge_with_maier) then begin
       good = where(raw.ha gt 0.0 and raw.ha_ew gt 0.0)
       data[good].h_alpha_continuum[0] = raw[good].ha/raw[good].ha_ew*1D-17
       data[good].h_alpha_continuum[1] = raw[good].ha_err*1D-17
       data[good].nii_6548_continuum  = data[good].h_alpha_continuum
       data[good].nii_6584_continuum  = data[good].h_alpha_continuum
    endif
       
; ---------------------------------------------------------------------------    
; EW's
; ---------------------------------------------------------------------------    

    good = where(raw.oii_ew gt 0.0)
    data[good].oii_3727_ew[0] = raw[good].oii_ew
    data[good].oii_3727_ew[1] = raw[good].oii_ew * raw[good].oii_err / raw[good].oii

    good = where(raw.hb_ew gt 0.0)
    data[good].h_beta_ew[0] = raw[good].hb_ew
    data[good].h_beta_ew[1] = raw[good].hb_ew * raw[good].hb_err / raw[good].hb

    good = where((data.oiii_5007[1] gt 0.0) and (data.oiii_5007_continuum[1] gt 0.0))
    data[good].oiii_5007_ew[0] = data[good].oiii_5007[0] / data[good].oiii_5007_continuum[0]
    data[good].oiii_5007_ew[1] = data[good].oiii_5007_ew[0] * raw[good].oiii_err / raw[good].oiii

    data[good].oiii_4959_ew = data[good].oiii_5007 / 3.0 ; assume this is true

    if keyword_set(merge_with_maier) then begin
       good = where(raw.ha_ew gt 0.0)
       data[good].h_alpha_ew[0] = raw[good].ha_ew
       data[good].h_alpha_ew[1] = raw[good].ha_ew * raw[good].ha_err / raw[good].ha

       good = where((data.nii_6584[1] gt 0.0) and (data.nii_6584_continuum[1] gt 0.0))
       data[good].nii_6584_ew[0] = data[good].nii_6584[0] / data[good].nii_6584_continuum[0]
       data[good].nii_6584_ew[1] = data[good].nii_6584_ew[0] * raw[good].nii_err / raw[good].nii

       data[good].nii_6548_ew = data[good].nii_6584 / 3.0 ; assume this is true
    endif
       
; ---------------------------------------------------------------------------    
; assign R23 branches
; ---------------------------------------------------------------------------    

    splog, 'Assigning R23 branches'
    data.r23branch = raw.branch
    
; ---------------------------------------------------------------------------    
; compute the reddening, metallicity, and other miscellaneous
; properties, then write out 
; ---------------------------------------------------------------------------    
    
;   idata = iunred_linedust(data,snrcut=snrcut,/silent,/nopropagate)
    ages_mz_log12oh, data, nmonte=500L, snrcut=0.0, final_ohdust=log12oh, $
      logfile=path+'mz_log12oh.03lilly.log'
    outdata = struct_addtags(data,log12oh)

    splog, 'Writing '+path+root+'.fits'
    mwrfits, outdata, path+root+'.fits', /create

;; just use the NII/Ha<-1.0 cut; there is no reason why Maier's
;; abundances are correct (plus look at his uncertainties!)
;    
;    adata = im_abundance(data,snrcut=snrcut)
;    outdata = struct_addtags(data,adata)
;    outdata = abundance_catalogs_log12oh(outdata);,ewbranch=raw.branch,branch=raw.branch)
;;   niceprint, outdata.galaxy, outdata.r23branch_kk04_ew, outdata.zstrong_ew_r23

; ---------------------------------------------------------------------------    
; write out statistics
; ---------------------------------------------------------------------------    

    indx = where((outdata.m_b gt -900.0) and (outdata.zstrong_ew_alpha_unity_12oh_kk04 gt -900.0),nindx)
    
    print, '##################################################'
    splog, 'Lilly et al. (2003): '+string(nindx,format='(I0.0)')+' galaxies.'
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
       
       im_window, 0, xratio=0.4, /square

       plot, outdata.m_b, outdata.zstrong_ew_alpha_unity_12oh_kk04, $
         xrange=[-23,-15], yrange=[8.0,9.2], ps=4, xsty=3, ysty=3, sym=2, $
         charsize=2.0, charthick=2.0
       cc = get_kbrd(1)

       niceprint, outdata.galaxy, outdata.mass_glaze, outdata.mass
       plot, outdata.mass_glaze, outdata.mass, xrange=[8,12.0], yrange=[8,12.0], $
         ps=4, xsty=3, ysty=3, sym=2, charsize=2.0, charthick=2.0
       oplot, !x.crange, !y.crange, thick=2.0
       jj = im_stats(outdata.mass_glaze-outdata.mass,/verbose)
       cc = get_kbrd(1)
       
       niceprint, outdata.galaxy, outdata.m_b_lilly, outdata.m_b
       plot, outdata.m_b_lilly, outdata.m_b, xrange=[-18.5,-23.0], $
         yrange=[-18.5,-23.0], ps=4, xsty=3, ysty=3, sym=2, charsize=2.0, charthick=2.0
       oplot, !x.crange, !y.crange, thick=2.0
       jj = im_stats(outdata.m_b_lilly-outdata.m_b,/verbose)
       cc = get_kbrd(1)

; compare the BRODWIN and SAVAGLIO photometry

       phot1 = rsex('cfrs_phot_from_savaglio.dat')
       phot2 = rsex('cfrs_phot_from_brodwin.dat')
       
       match, phot1.galaxy, phot2.galaxy, indx1, indx2
;   niceprint, phot1[indx1].galaxy, phot2[indx2].galaxy

       plot, phot1[indx1].i_iso, phot2[indx2].i, ps=4, xr=[21,23], yr=[21,23], $
         sym=2, charsize=2.0, charthick=2.0
       oplot, !x.crange, !y.crange, thick=2
       cc = get_kbrd(1)
       
       plot, phot1[indx1].i_iso+(phot1[indx1].b3-phot1[indx1].i3), $
         phot2[indx2].i+(phot2[indx2].b_aper-phot2[indx2].i_aper), ps=4, $
         xr=[21.5,25], yr=[21.5,25], sym=2, charsize=2.0, charthick=2.0
       oplot, !x.crange, !y.crange, thick=2
       cc = get_kbrd(1)
       
       plot, phot1[indx1].i_iso+(phot1[indx1].v3-phot1[indx1].i3), $
         phot2[indx2].i+(phot2[indx2].v_aper-phot2[indx2].i_aper), ps=4, $
         xr=[21.5,25], yr=[21.5,25], sym=2, charsize=2.0, charthick=2.0
       oplot, !x.crange, !y.crange, thick=2
       cc = get_kbrd(1)
       
       plot, phot1[indx1].i_iso+(phot1[indx1].k3-phot1[indx1].i3), $
         phot2[indx2].i+(phot2[indx2].k_aper-phot2[indx2].i_aper), ps=4, $
         xr=[17,24], yr=[17,24], sym=2, charsize=2.0, charthick=2.0
       oplot, !x.crange, !y.crange, thick=2
       cc = get_kbrd(1)
    
    endif
       
return
end    
