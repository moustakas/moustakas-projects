pro lcs_xray

    cl = rsex(getenv('IDL_PROJECTS_DIR')+'/lcs/lcs_sample.cat')
    xx = mrdfits(getenv('CATALOGS_DIR')+'/mcxc_xray_catalog/mcxc.fits.gz',1)

    spherematch, xx.raj2000, xx.dej2000, cl.ra, $
      cl.dec, 360D/3600.0, m1, m2    

    lcspath = getenv('LCS_DATA')
    im_mwrfits, xx[m1], lcspath+'lcs_mcxc.fits', /clobber
    
return
end
    
