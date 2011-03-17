pro ages_match_twomass
; jm09may15nyu - match the AGES and 2MASS catalogs

    common match_twomass, twomass

    if (n_elements(twomass) eq 0) then twomass = $
      hogg_mrdfits(vagc_name('object_twomass'),1,nrow=28800L,$
      columns=['ra','decl','twomass_tag','k_m_ext','k_msig_ext',$
      'h_m_ext','h_msig_ext','j_m_ext','j_msig_ext'])

    ages = mrdfits(ages_path(/catalogs)+'catalog.cat.noguidestars.fits.gz',1)
    ages.ra = ages.ra*15.0D

    spherematch, twomass.ra, twomass.decl, ages.ra, $
      ages.dec, 1.5/3600.0, m1, m2
    srt = sort(m2)
    m1 = m1[srt]
    m2 = m2[srt]

    tmass = im_empty_structure(twomass[0],empty_value=-999.0,$
      ncopies=n_elements(ages))
    tmass[m2] = twomass[m1]
        
    im_mwrfits, tmass, ages_path(/mycatalogs)+'ages.twomass.phot.fits'

return
end
    
