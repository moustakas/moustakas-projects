pro write_integrated_abundances_sample, atlasdust, atlasnodust, nfgsdust, nfgsnodust, write=write, _extra=extra
; jm05aug26uofa - write integrated spectroscopic sample for the
;                 AGES/MZ paper 

    outpath = atlas_path(/projects)+'abundances/'
    snrcut = 3.0
    
    if (n_elements(atlasdust) eq 0L) then atlasdust = read_integrated(atlasnodust=atlasnodust)
    if (n_elements(nfgsdust) eq 0L) then nfgsdust = read_nfgs(nfgsnodust=nfgsnodust)

; append the two samples    
    
    intdust = struct_append(atlasdust,nfgsdust)
    intnodust = struct_append(atlasnodust,nfgsnodust)
    ngalaxy = n_elements(intdust)

; define more structure tags    
    
    abund = where((intdust.m_b_obs gt -900.0) and $
      (intdust.h_beta[0]/intdust.h_beta[1] gt snrcut) and $
      (intdust.h_alpha[0]/intdust.h_alpha[1] gt snrcut) and $
      (intdust.nii_6584[0]/intdust.nii_6584[1] gt snrcut) and $
      (intdust.oii_3727[0]/intdust.oii_3727[1] gt snrcut) and $
      (intdust.oiii_5007[0]/intdust.oiii_5007[1] gt snrcut),nabund)

    hii = where(strtrim(intdust[abund].bpt_pure_nii_class,2) eq 'HII',nhii)
    agn = where(strtrim(intdust[abund].bpt_pure_nii_class,2) eq 'AGN',nagn)
    splog, 'Abundances sample: '+string(nabund,format='(I0)')+'/'+string(ngalaxy,format='(I0)')+' ('+$
      string(round(100.0*nabund/float(ngalaxy)),format='(I0)')+'%).'
    splog, 'Abundances sample (star-forming): '+string(nhii,format='(I0)')+'/'+string(nabund,format='(I0)')+' ('+$
      string(round(100.0*nhii/float(nabund)),format='(I0)')+'%).'
    splog, 'Abundances sample (type 2 AGN): '+string(nagn,format='(I0)')+'/'+string(nabund,format='(I0)')+' ('+$
      string(round(100.0*nagn/float(nabund)),format='(I0)')+'%).'
    print

; assign galaxies to the appropriate R23 branch

    splog, 'Assigning R23 branches.'
    ohnodust_kk04 = im_assign_r23branch(intnodust,branchmethod=2L,mb=intdust.m_b_obs,subindx=abund[hii],/kk04,_extra=extra)
    ohnodust_m91  = im_assign_r23branch(intnodust,branchmethod=2L,mb=intdust.m_b_obs,subindx=abund[hii],/m91,_extra=extra)
    ohnodust_pt05 = im_assign_r23branch(intnodust,branchmethod=2L,mb=intdust.m_b_obs,subindx=abund[hii],/pt05,_extra=extra)
    ohnodust = struct_addtags(struct_addtags(ohnodust_pt05,ohnodust_m91),ohnodust_kk04)
    print
    
; abundance statistics: corrected fluxes

    pt05_good = where((ohnodust[abund[hii]].zstrong_12oh_pt05 gt -900.0),npt05_good)
    pt05_good_upper = where((strmatch(ohnodust[abund[hii[pt05_good]]].r23branch_pt05,'U_*')),npt05_good_upper)
    pt05_good_lower = where((strmatch(ohnodust[abund[hii[pt05_good]]].r23branch_pt05,'L_*')),npt05_good_lower)
    pt05_good_ambig = where((strmatch(ohnodust[abund[hii[pt05_good]]].r23branch_pt05,'A')),npt05_good_ambig)
    splog, 'PT05 abundances (corrected): '+string(npt05_good,format='(I0)')+'/'+string(nhii,format='(I0)')+' ('+$
      string(round(100.0*npt05_good/nhii),format='(I0)')+'%).'
    splog, '   Upper: '+string(npt05_good_upper,format='(I0)')+'/'+string(npt05_good,format='(I0)')+' ('+$
      string(round(100.0*npt05_good_upper/npt05_good),format='(I0)')+'%).'
    splog, '   Lower: '+string(npt05_good_lower,format='(I0)')+'/'+string(npt05_good,format='(I0)')+' ('+$
      string(round(100.0*npt05_good_lower/npt05_good),format='(I0)')+'%).'
    splog, '   Ambig: '+string(npt05_good_ambig,format='(I0)')+'/'+string(npt05_good,format='(I0)')+' ('+$
      string(round(100.0*npt05_good_ambig/npt05_good),format='(I0)')+'%).'
    print

    m91_good = where((ohnodust[abund[hii]].zstrong_12oh_m91 gt -900.0),nm91_good)
    m91_good_upper = where((strmatch(ohnodust[abund[hii[m91_good]]].r23branch_m91,'U_*')),nm91_good_upper)
    m91_good_lower = where((strmatch(ohnodust[abund[hii[m91_good]]].r23branch_m91,'L_*')),nm91_good_lower)
    m91_good_ambig = where((strmatch(ohnodust[abund[hii[m91_good]]].r23branch_m91,'A')),nm91_good_ambig)
    splog, 'M91 abundances (corrected): '+string(nm91_good,format='(I0)')+'/'+string(nhii,format='(I0)')+' ('+$
      string(round(100.0*nm91_good/nhii),format='(I0)')+'%).'
    splog, '   Upper: '+string(nm91_good_upper,format='(I0)')+'/'+string(nm91_good,format='(I0)')+' ('+$
      string(round(100.0*nm91_good_upper/nm91_good),format='(I0)')+'%).'
    splog, '   Lower: '+string(nm91_good_lower,format='(I0)')+'/'+string(nm91_good,format='(I0)')+' ('+$
      string(round(100.0*nm91_good_lower/nm91_good),format='(I0)')+'%).'
    splog, '   Ambig: '+string(nm91_good_ambig,format='(I0)')+'/'+string(nm91_good,format='(I0)')+' ('+$
      string(round(100.0*nm91_good_ambig/nm91_good),format='(I0)')+'%).'
    print

    kk04_good = where((ohnodust[abund[hii]].zstrong_12oh_kk04 gt -900.0),nkk04_good)
    kk04_good_upper = where((strmatch(ohnodust[abund[hii[kk04_good]]].r23branch_kk04,'U_*')),nkk04_good_upper)
    kk04_good_lower = where((strmatch(ohnodust[abund[hii[kk04_good]]].r23branch_kk04,'L_*')),nkk04_good_lower)
    kk04_good_ambig = where((strmatch(ohnodust[abund[hii[kk04_good]]].r23branch_kk04,'A')),nkk04_good_ambig)
    splog, 'KK04 abundances (corrected): '+string(nkk04_good,format='(I0)')+'/'+string(nhii,format='(I0)')+' ('+$
      string(round(100.0*nkk04_good/nhii),format='(I0)')+'%).'
    splog, '   Upper: '+string(nkk04_good_upper,format='(I0)')+'/'+string(nkk04_good,format='(I0)')+' ('+$
      string(round(100.0*nkk04_good_upper/nkk04_good),format='(I0)')+'%).'
    splog, '   Lower: '+string(nkk04_good_lower,format='(I0)')+'/'+string(nkk04_good,format='(I0)')+' ('+$
      string(round(100.0*nkk04_good_lower/nkk04_good),format='(I0)')+'%).'
    splog, '   Ambig: '+string(nkk04_good_ambig,format='(I0)')+'/'+string(nkk04_good,format='(I0)')+' ('+$
      string(round(100.0*nkk04_good_ambig/nkk04_good),format='(I0)')+'%).'
    print, '--------------------------------------------------'
    print

    if keyword_set(write) then begin

       splog, 'Appending additional tags.'
;      intdust = struct_addtags(temporary(intdust),struct_addtags(moretags,ohdust))
       intnodust = struct_addtags(temporary(intnodust),struct_addtags(moretags,ohnodust))

; --------------------
       
;      splog, 'Writing '+outpath+'integrated_abundances_speclinefit.fits.gz'
;      mwrfits, intdust[abund], outpath+'integrated_abundances_speclinefit.fits', /create
;      spawn, ['gzip -f '+outpath+'integrated_abundances_speclinefit.fits'], /sh
;
;      splog, 'Writing '+outpath+'integrated_abundances_speclinefit_nodust.fits.gz'
;      mwrfits, intnodust[abund], outpath+'integrated_abundances_speclinefit_nodust.fits', /create
;      spawn, ['gzip -f '+outpath+'integrated_abundances_speclinefit_nodust.fits'], /sh

; --------------------
       
       splog, 'Writing '+outpath+'integrated_abundances_hii_speclinefit.fits.gz'
       mwrfits, intdust[abund[hii]], outpath+'integrated_abundances_hii_speclinefit.fits', /create
       spawn, ['gzip -f '+outpath+'integrated_abundances_hii_speclinefit.fits'], /sh

       splog, 'Writing '+outpath+'integrated_abundances_hii_speclinefit_nodust.fits.gz'
       mwrfits, intnodust[abund[hii]], outpath+'integrated_abundances_hii_speclinefit_nodust.fits', /create
       spawn, ['gzip -f '+outpath+'integrated_abundances_hii_speclinefit_nodust.fits'], /sh

; --------------------

;      splog, 'Writing '+outpath+'integrated_abundances_agn_speclinefit.fits.gz'
;      mwrfits, intdust[abund[agn]], outpath+'integrated_abundances_agn_speclinefit.fits', /create
;      spawn, ['gzip -f '+outpath+'integrated_abundances_agn_speclinefit.fits'], /sh
;
;      splog, 'Writing '+outpath+'integrated_abundances_agn_speclinefit_nodust.fits.gz'
;      mwrfits, intnodust[abund[agn]], outpath+'integrated_abundances_agn_speclinefit_nodust.fits', /create
;      spawn, ['gzip -f '+outpath+'integrated_abundances_agn_speclinefit_nodust.fits'], /sh

    endif

return
end
    
