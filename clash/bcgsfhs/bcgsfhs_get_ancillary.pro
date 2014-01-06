pro bcgsfhs_get_ancillary

    path = bcgsfhs_path(/ancillary)
    
    dr9 = mrdfits('photoPosPlate-dr9.fits.gz',1)
    
    ss = read_bcgsfhs_sample()
    spherematch, dr9.ra, dr9.dec, 15D*hms2dec(ss.ra), hms2dec(ss.dec), 30D/3600.0, m1, m2
    galex = mrdfits('sdss_dr9_galex_gr6.fits.gz',1,rows=m1)

    out = struct_addtags(struct_trimtags(ss[m2],select=['cluster','shortname']),dr9[m1])
    im_mwrfits, out, path+'ancillary_sdss_dr9.fits', /clob

    out = struct_addtags(struct_trimtags(ss[m2],select=['cluster','shortname']),galex)
    im_mwrfits, out, path+'ancillary_galex_gr6.fits', /clob
    
    wise = mrdfits('sdss_dr9_wise.fits.gz',1,rows=m1)
    im_mwrfits, out, path+'ancillary_wise.fits', /clob

return
end
    
