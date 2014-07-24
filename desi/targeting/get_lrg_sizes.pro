pro get_lrg_sizes

    path = getenv('IM_PROJECTS_DIR')+'/desi/lrgs/'
    readcol, path+'lrgcoords.txt', ra, dec, format='D,D'
    ngal = n_elements(ra)

; read the ACS-GC catalog from Griffith+13
    acs = mrdfits(getenv('IM_DATA_DIR')+'/acs-gc/cosmos_i_public_catalog_V1.0.fits.gz',1)
    keep = where(acs.flag_galfit_hi eq 0)
    acs = acs[keep]

    out = replicate({ra: 0D, dec: 0D, acs_ra: -999D, acs_dec: -999D, $
      objno: -999L, re: -999.0, n: -999.0},ngal)
    out.ra = ra
    out.dec = dec

    spherematch, out.ra, out.dec, acs.ra, acs.dec, 1D/3600.0, m1, m2
    out[m1].acs_ra = acs[m2].ra
    out[m1].acs_dec = acs[m2].dec
    out[m1].re = acs[m2].re_galfit_hi*0.05 ; [arcsec]
    out[m1].n = acs[m2].n_galfit_hi ; [arcsec]
    out[m1].objno = acs[m2].objno

    out_acs = im_empty_structure(acs[0],ncopies=ngal)
    out_acs[m1] = acs[m2]
    
    need = where(out.re lt 0,nneed)
    spherematch, out[need].ra, out[need].dec, acs.ra, acs.dec, 1D/3600.0, m1, m2
    if m1[0] ne -1 then begin
       out[need[m1]].re = acs[m2].re_galfit_lo*0.05 ; [arcsec]
       out[need[m1]].n = acs[m2].n_galfit_lo*0.05   ; [arcsec]
    endif

    im_mwrfits, out, path+'lrg_sizes.fits', /clobber
    im_mwrfits, out_acs, path+'lrg_sizes_acs.fits', /clobber
    struct_print, out
    
stop    

return
end
    
