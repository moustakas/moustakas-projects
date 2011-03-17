;+
; NAME:
;   BUILD_SG1120_PARENT_CATALOG
;
; PURPOSE:
;   Merge all the SG1120 photometric catalogs together.
;
; ToDo:
;   * Add HST, FLAMINGOS.
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2009 Jun 15, NYU - based on older code
;-

pro build_sg1120_parent_catalog, phot

    common sg1120_phot, ebv1, hst1, bcat1, vcat1, $
      rcat1, gprimecat1, rprimecat1, kscat1, sdss1
    
; redshift catalog; the suffix specifies which 
    zcat1 = sg1120_read_zcat(suffix=suffix)
    
    if (n_elements(sdss1) eq 0L) then begin
       sexpath = sg1120_path(/sex)
       splog, 'Reading '+sexpath+'sg1120_sdss_dr7.fits.gz'
       sdss1 = mrdfits(sexpath+'sg1120_sdss_dr7.fits.gz',1)
    endif

; VIMOS catalogs
    vversion = vimos_catalogs_version()
    vcatpath = vimos_path(/catalogs)
    vcatnames = vcatpath+'sg1120_'+['B','V','R']+'_'+vversion+'.chi2.cat'
    if (n_elements(bcat1) eq 0L) then begin
       splog, 'Reading '+vcatnames[0]
       bcat1 = rsex(vcatnames[0])
    endif
    if (n_elements(vcat1) eq 0L) then begin
       splog, 'Reading '+vcatnames[1]
       vcat1 = rsex(vcatnames[1])
    endif
    if (n_elements(rcat1) eq 0L) then begin
       splog, 'Reading '+vcatnames[2]
       rcat1 = rsex(vcatnames[2])
    endif
    
; LDSS3 catalogs
    lversion = ldss3_catalogs_version()
    lcatpath = ldss3_path(/catalogs)
    lcatnames = lcatpath+'sg1120_'+['gprime','rprime']+'_'+lversion+'.chi2.cat'
    if (n_elements(gprimecat1) eq 0L) then begin
       splog, 'Reading '+lcatnames[0]
       gprimecat1 = rsex(lcatnames[0])
    endif
    if (n_elements(rprimecat1) eq 0L) then begin
       splog, 'Reading '+lcatnames[1]
       rprimecat1 = rsex(lcatnames[1])
    endif

; Ks catalog
    kversion = flamingos_catalogs_version()
    kcatpath = flamingos_path(/catalogs)
    kcatnames = kcatpath+'sg1120_Ks_'+kversion+'.cat'
    if (n_elements(kscat1) eq 0L) then begin
       splog, 'Reading '+kcatnames[0]
       kscat1 = rsex(kcatnames[0])
    endif

; HST catalog
    hcatpath = hst_path(/catalogs)
    hcatname = hcatpath+'sg1120_hst.cat'
    if (n_elements(hst1) eq 0L) then begin
       splog, 'Reading '+hcatname
       hst1 = rsex(hcatname)
    endif

; initialize the output data structure
    phot = {$
      id:                  0L, $ ; unique ID number, zero-indexed
      ra:                0.0D, $
      dec:               0.0D, $
      runid:               '', $
      z:                  0.0, $
      q:                    0, $
      ebv_mw:          -999.0, $
      b_id:             -999L, $
      v_id:             -999L, $
      r_id:             -999L, $
      gprime_id:        -999L, $
      rprime_id:        -999L, $
      ks_id:           -999L, $
      f814_id:          -999L, $
      sdss_type:           -1, $
      bvr_flags:      [0,0,0], $
      bvr_flag_split: [0,0,0], $
      gr_flags:         [0,0], $
      gr_flag_split:    [0,0], $
      ks_flags:             0, $
      ks_flag_split:        0, $
      f814_flags:           0, $
      f814_pointing:        0, $
      phot_b:          -999.0, $
      phot_b_err:      -999.0, $
      phot_v:          -999.0, $
      phot_v_err:      -999.0, $
      phot_r:          -999.0, $
      phot_r_err:      -999.0, $
      phot_gprime:     -999.0, $
      phot_gprime_err: -999.0, $
      phot_rprime:     -999.0, $
      phot_rprime_err: -999.0, $
      phot_sdssu:      -999.0, $
      phot_sdssu_err:  -999.0, $
      phot_sdssg:      -999.0, $
      phot_sdssg_err:  -999.0, $
      phot_sdssr:      -999.0, $
      phot_sdssr_err:  -999.0, $
      phot_sdssi:      -999.0, $
      phot_sdssi_err:  -999.0, $
      phot_ks:         -999.0, $
      phot_ks_err:     -999.0, $
      phot_f814:       -999.0, $
      phot_f814_err:   -999.0}

; the VIMOS R-band catalog defines the parent source list
    ngal = n_elements(rcat1)
    phot = replicate(phot,ngal)
    
    phot.id = lindgen(ngal)
    phot.ra = rcat1.xwin_world
    phot.dec = rcat1.ywin_world

    if (n_elements(ebv1) eq 0L) then begin
       splog, 'Computing MW E(B-V)'
       glactc, phot.ra, phot.dec, 2000.0, gl, gb, 1, /degree
       ebv1 = dust_getval(gl,gb,/interp)
    endif
    phot.ebv_mw = ebv1

; VIMOS; even though the VIMOS/BVR catalogs are defined using the same
; detection image, spherematch anyway to ensure consistency
    splog, 'Sphere-matching VIMOS/B'
    spherematch, phot.ra, phot.dec, bcat1.alpha_j2000, $
      bcat1.delta_j2000, 1.0/3600.0, m1, m2
    phot[m1].b_id = bcat1[m2].number
    phot[m1].phot_b = bcat1[m2].mag_auto
    phot[m1].phot_b_err = bcat1[m2].magerr_auto
    phot[m1].bvr_flags[0] = bcat1[m2].flags
    phot[m1].bvr_flag_split[0] = bcat1[m2].flag_split
    
    splog, 'Sphere-matching VIMOS/V'
    spherematch, phot.ra, phot.dec, vcat1.alpha_j2000, $
      vcat1.delta_j2000, 1.0/3600.0, m1, m2
    phot[m1].v_id = vcat1[m2].number
    phot[m1].phot_v = vcat1[m2].mag_auto
    phot[m1].phot_v_err = vcat1[m2].magerr_auto
    phot[m1].bvr_flags[1] = vcat1[m2].flags
    phot[m1].bvr_flag_split[1] = vcat1[m2].flag_split
    
    splog, 'Sphere-matching VIMOS/R'
    spherematch, phot.ra, phot.dec, rcat1.alpha_j2000, $
      rcat1.delta_j2000, 1.0/3600.0, m1, m2
    phot[m1].r_id = rcat1[m2].number
    phot[m1].phot_r = rcat1[m2].mag_auto
    phot[m1].phot_r_err = rcat1[m2].magerr_auto
    phot[m1].bvr_flags[2] = rcat1[m2].flags
    phot[m1].bvr_flag_split[2] = rcat1[m2].flag_split

; LDSS3
    splog, 'Sphere-matching LDSS3/g'
    spherematch, phot.ra, phot.dec, gprimecat1.alpha_j2000, $
      gprimecat1.delta_j2000, 1.0/3600.0, m1, m2
    phot[m1].gprime_id = gprimecat1[m2].number
    phot[m1].phot_gprime = gprimecat1[m2].mag_auto
    phot[m1].phot_gprime_err = gprimecat1[m2].magerr_auto
    phot[m1].gr_flags[0] = gprimecat1[m2].flags
    phot[m1].gr_flag_split[0] = gprimecat1[m2].flag_split
    
    splog, 'Sphere-matching LDSS3/r'
    spherematch, phot.ra, phot.dec, rprimecat1.alpha_j2000, $
      rprimecat1.delta_j2000, 1.0/3600.0, m1, m2
    phot[m1].rprime_id = rprimecat1[m2].number
    phot[m1].phot_rprime = rprimecat1[m2].mag_auto
    phot[m1].phot_rprime_err = rprimecat1[m2].magerr_auto
    phot[m1].gr_flags[1] = rprimecat1[m2].flags
    phot[m1].gr_flag_split[1] = rprimecat1[m2].flag_split

; FLAMINGOS
    splog, 'Sphere-matching FLAMINGOS/Ks'
    spherematch, phot.ra, phot.dec, kscat1.alpha_j2000, $
      kscat1.delta_j2000, 1.0/3600.0, m1, m2
    phot[m1].ks_id = kscat1[m2].number
    phot[m1].phot_ks = kscat1[m2].mag_auto
    phot[m1].phot_ks_err = kscat1[m2].magerr_auto
    phot[m1].ks_flags = kscat1[m2].flags
    phot[m1].ks_flag_split = kscat1[m2].flag_split
    
; HST    
    splog, 'Sphere-matching HST/F814'
    spherematch, phot.ra, phot.dec, hst1.alpha_j2000, $
      hst1.delta_j2000, 1.0/3600.0, m1, m2
    phot[m1].f814_id = hst1[m2].number
    phot[m1].phot_f814 = hst1[m2].mag_auto
    phot[m1].phot_f814_err = hst1[m2].magerr_auto
    phot[m1].f814_flags = hst1[m2].flags
    phot[m1].f814_pointing = hst1[m2].pointing
    
; set crummy photometry to -999.0
    crap = where((phot.phot_b le 0.0) or (phot.phot_b ge 90.0) or (phot.phot_b_err ge 1.0),ncrap)
    if (ncrap ne 0) then begin
       phot[crap].phot_b = -999.0
       phot[crap].phot_b_err = -999.0
    endif
    crap = where((phot.phot_v le 0.0) or (phot.phot_v ge 90.0) or (phot.phot_v_err ge 1.0),ncrap)
    if (ncrap ne 0) then begin
       phot[crap].phot_v = -999.0
       phot[crap].phot_v_err = -999.0
    endif
    crap = where((phot.phot_r le 0.0) or (phot.phot_r ge 90.0) or (phot.phot_r_err ge 1.0),ncrap)
    if (ncrap ne 0) then begin
       phot[crap].phot_r = -999.0
       phot[crap].phot_r_err = -999.0
    endif
    crap = where((phot.phot_gprime le 0.0) or (phot.phot_gprime ge 90.0) or (phot.phot_gprime_err ge 1.0),ncrap)
    if (ncrap ne 0) then begin
       phot[crap].phot_gprime = -999.0
       phot[crap].phot_gprime_err = -999.0
    endif
    crap = where((phot.phot_rprime le 0.0) or (phot.phot_rprime ge 90.0) or (phot.phot_rprime_err ge 1.0),ncrap)
    if (ncrap ne 0) then begin
       phot[crap].phot_rprime = -999.0
       phot[crap].phot_rprime_err = -999.0
    endif
    crap = where((phot.phot_f814 le 0.0) or (phot.phot_f814 ge 90.0) or (phot.phot_f814_err ge 1.0),ncrap)
    if (ncrap ne 0) then begin
       phot[crap].phot_f814 = -999.0
       phot[crap].phot_f814_err = -999.0
    endif

; SDSS
    splog, 'Sphere-matching against SDSS'
    spherematch, phot.ra, phot.dec, sdss1.ra, sdss1.dec, $
      1.0/3600.0, m1, m2
    phot[m1].sdss_type = sdss1[m2].type
    phot[m1].phot_sdssu = sdss1[m2].u
    phot[m1].phot_sdssu_err = sdss1[m2].uerr
    phot[m1].phot_sdssg = sdss1[m2].g
    phot[m1].phot_sdssg_err = sdss1[m2].gerr
    phot[m1].phot_sdssr = sdss1[m2].r
    phot[m1].phot_sdssr_err = sdss1[m2].rerr
    phot[m1].phot_sdssi = sdss1[m2].i
    phot[m1].phot_sdssi_err = sdss1[m2].ierr

; redshift catalog
    splog, 'Sphere-matching the redshift catalog'
    spherematch, phot.ra, phot.dec, zcat1.ra, zcat1.dec, $
      1.5/3600.0, m1, m2
    phot[m1].z = zcat1[m2].z
    phot[m1].q = zcat1[m2].q
    phot[m1].runid = zcat1[m2].runid
    
; write out
    outpath = sg1120_path(/analysis)
    version = sg1120_version(/parent)
    outfile = outpath+'sg1120_parent_'+suffix+'_'+version+'.fits'
;   outfile = outpath+'sg1120_parent_'+suffix+'.fits'
    im_mwrfits, phot, outfile

stop    
    
return
end    
