pro build_redmapper_photometry, out
; jm13mar28siena - merge the GALEX, SDSS, and WISE photometry
; jm14aug06siena - updated to v5.10, which no longer includes GALEX
;   photometry and now includes unWISE photometry 

; MEM_MATCH_ID: unique id
; RA, DEC: most likely central galaxy position
; LAMBDA_CHISQ(_E): richness (& error)
; Z_LAMBDA(_E): redmapper photo-z (& Gaussian error)
; BCG_SPEC_Z: spectroscopic redshift of the central galaxy if available
; PZBINS, PZ: P(z) for the cluster, at 21 points.  (peak is at Z_LAMBDA).
; 
; In the members list, they are:
; MEM_MATCH_ID: id keyed to cluster list
; RA, DEC: position
; P: probability of membership (P_mem in the paper)
; R: radius (h^-1 Mpc from center)

; In the members list, there is a ".Z" tag that is the same as the
; .Z_LAMBDA tag for the cluster.
    
    filt = redmapper_filterlist()
    nfilt = n_elements(filt)
    
    path = redmapper_path(version=ver)
;   bcgs = mrdfits(path+'dr8_run_redmapper_'+ver+$
;     '_lgt20_catalog.fits.gz',1)
    cat = mrdfits(path+'dr8_run_redmapper_'+ver+$
      '_lgt5_catalog_members.fits.gz',1);,range=[0,1000])
    ngal = n_elements(cat)

    sdss = mrdfits(path+'redmapper_'+ver+'_sdss.fits.gz',1);,range=[0,1000])
    unwise = mrdfits(path+'redmapper_'+ver+'_unwise.fits.gz',1);,range=[0,1000])
;   galex = mrdfits(path+'redmapper_'+ver+'_galex_gr6.fits.gz',1)

;   sdss_to_maggies, smaggies, sivarmaggies, cas=sdss, flux='cmodel'
    sdss_to_maggies, modelmaggies, modelivarmaggies, cas=sdss, flux='model'
    sdss_to_maggies, cmodelmaggies, cmodelivarmaggies, cas=sdss, flux='cmodel'
    ratio = cmodelmaggies[2,*]/modelmaggies[2,*]
    neg = where(modelmaggies[2,*] le 0)
;   if (neg[0] ne -1L) then message, 'Bad!'
    if (neg[0] ne -1L) then ratio[neg] = 1.0

    factor = rebin(ratio,5,ngal)
    smaggies = modelmaggies*factor
    sivarmaggies = modelivarmaggies/factor^2

    unwise_to_maggies, unwise, wmaggies, wivarmaggies
;   wise_to_maggies, wise, wmaggies, wivarmaggies, /mpro
;   im_galex_to_maggies, galex, gmaggies, givarmaggies

;   out = struct_trimtags(cat)
;   out = struct_trimtags(cat,select=['mem_match_id','ra','dec','z',$
;     'r','p','photoid'])
    out = struct_addtags(cat,replicate({$
      isbcg: 0,                   $
      maggies: fltarr(nfilt),     $
      ivarmaggies: fltarr(nfilt)},ngal))
;     irmaggies: fltarr(2),       $
;     irivarmaggies: fltarr(2)},ngal))

    out.maggies = [smaggies,wmaggies]
    out.ivarmaggies = [sivarmaggies,wivarmaggies]
;   out.maggies = [gmaggies,smaggies,wmaggies[0:1,*]]
;   out.ivarmaggies = [givarmaggies,sivarmaggies,wivarmaggies[0:1,*]]
;   out.irmaggies = wmaggies[2:3,*]
;   out.irivarmaggies = wivarmaggies[2:3,*]

; identify the centrals/BCGs
    out[where(cat.r eq 0)].isbcg = 1

; write out    
    im_mwrfits, out, path+'redmapper_'+ver+'_phot.fits', /clobber
    
return
end
    
