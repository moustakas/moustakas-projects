pro build_redmapper_photometry, out
; jm13mar28siena - merge the GALEX, SDSS, and WISE photometry 

    filt = redmapper_filterlist()
    nfilt = n_elements(filt)
    
    path = redmapper_path(/catalogs,version=ver)
;   bcgs = mrdfits(path+'dr8_run_redmapper_'+ver+$
;     '_lgt20_catalog.fits.gz',1)
    cat = mrdfits(path+'dr8_run_redmapper_'+ver+$
      '_lgt20_catalog_members.fits.gz',1)
    ngal = n_elements(cat)

    sdss = mrdfits(path+'redmapper_'+ver+'_sdss.fits.gz',1)
    galex = mrdfits(path+'redmapper_'+ver+'_galex_gr6.fits.gz',1)
    wise = mrdfits(path+'redmapper_'+ver+'_wise.fits.gz',1)

    wise_to_maggies, wise, wmaggies, wivarmaggies, /mpro
    im_galex_to_maggies, galex, gmaggies, givarmaggies
    sdss_to_maggies, smaggies, sivarmaggies, cas=sdss, flux='cmodel'

    out = struct_trimtags(cat,select=['mem_match_id','ra','dec','z',$
      'r','p','photoid','isbcg'])
    out = struct_addtags(out,replicate({isbcg: 0, maggies: fltarr(nfilt), $
      ivarmaggies: fltarr(nfilt), irmaggies: fltarr(2), $
      irivarmaggies: fltarr(2)},ngal))
    out.maggies = [gmaggies,smaggies,wmaggies[0:1,*]]
    out.ivarmaggies = [givarmaggies,sivarmaggies,wivarmaggies[0:1,*]]
    out.irmaggies = wmaggies[2:3,*]
    out.irivarmaggies = wivarmaggies[2:3,*]

; identify the centrals/BCGs
    out[where(cat.r eq 0)].isbcg = 1

; write out    
    im_mwrfits, out, path+'redmapper_'+ver+'_phot.fits', /clobber
    
return
end
    
