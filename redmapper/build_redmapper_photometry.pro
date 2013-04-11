pro build_redmapper_photometry, out
; jm13mar28siena - merge the GALEX, SDSS, and WISE photometry 

    path = redmapper_path(/catalogs,version=ver)

    filt = redmapper_filterlist()
    nfilt = n_elements(filt)
    
;   bcgs = mrdfits(path+'dr8_run_redmapper_'+ver+$
;     '_lgt20_catalog.fits.gz',1)
    cat = mrdfits(path+'dr8_run_redmapper_'+ver+$
      '_lgt20_catalog_members.fits.gz',1)
    ngal = n_elements(cat)

    sdss = mrdfits(path+'redmapper_v5.2_sdss.fits.gz',1)
    galex = mrdfits(path+'redmapper_v5.2_galex_gr6.fits.gz',1)
    wise = mrdfits(path+'redmapper_v5.2_wise.fits.gz',1)

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

; this misses all but four of the BCGs    
    splog, 'Finding BCGs...'
    out[where(cat.r eq 0)].isbcg = 1
    
;; do the unique ones first
;    all = lindgen(ngal)
;    uu = uniq(cat.photoid,sort(cat.photoid))
;
;    nonuu = all
;    remove, uu, nonuu
;    
;;   for ii = 0L, n_elements(bcgs)-1L do out[where(cat.photoid eq bcgs[ii].photoid)].isbcg = 1
;    these = cmset_op(cat.photoid,'and',bcgs.photoid,/index)
;    help, where(bcgs.photoid eq cat[these[0]].photoid)             
;    out[these].isbcg = 1

    im_mwrfits, out, path+'redmapper_'+ver+'_photometry.fits', /clobber
    
return
end
    
