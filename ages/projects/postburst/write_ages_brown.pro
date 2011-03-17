pro write_ages_brown
; jm07nov16nyu - write out an AGES data structure for the AGES
;                postburst/AGN project with M. Brown

    anc = read_ages(specdata=spec,/ispec)

    anc_tags = ['GALAXY','AGES_ID','RA','DEC','PASS','APER','Z','SPEC_WEIGHT',$;'MAIN_FLAG',$
      'X_FLUX','X_FLUX_ERR','X_LUM','X_HR','X_HR_ERR','KCORR_MASS','KCORR_CHI2',$
      'M_'+['U','B','V','R'],'M_'+['U','B','V','R']+'_IVAR']
    spec_tags = ['H_ALPHA_EW','H_BETA_EW','OII_3727_EW','NEV_3426','OIII_5007','NII_6584','H_ALPHA','H_BETA',$
      'D4000_NARROW','LICK_HD_A','LICK_HG_A',['D4000_NARROW','LICK_HD_A','LICK_HG_A']+'_MODEL']
      
    out = struct_addtags(struct_trimtags(anc,select=anc_tags),$
      struct_trimtags(spec,select=spec_tags))

    keep = where((anc.z gt 0.01) and (anc.z lt 1.0) and $
      (anc.fluxing_problem eq 0B) and (anc.main_flag eq 1B))
    outpath = ages_path(/projects)+'postburst/brown/'
    mwrfits, out[keep], outpath+'ages_brown.fits', /create
    spawn, 'gzip -f '+outpath+'ages_brown.fits', /sh

stop    
    
return
end
    
