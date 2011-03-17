pro oiinestats
; jm09apr01nyu - compute and write out the distribution of electron
;   densities in each redshift bin; no plots are generated

    oiinepath = deep2_path(/projects)+'oiine/'
    idlpath = getenv('HOME')+'/idl/projects/deep2/projects/oiine/'

    deepkcorr = read_oiine_sample(/ancillary)
    deepispec = read_oiine_sample(/ispec)
    sdsskcorr = read_oiine_sample(/ancillary,/sdss)
    sdssispec = read_oiine_sample(/ispec,/sdss)

; initialize the output data structure
    stats = {$
      z:           0.0,$
      zerr:        0.0,$
      ne_mean:     0.0,$
      ne_p05:      0.0,$
      ne_p16:      0.0,$
      ne_p25:      0.0,$
      ne_p50:      0.0,$
      ne_p84:      0.0,$
      ne_p75:      0.0,$
      ne_p95:      0.0,$
      ngal:         0L}
    quant = [0.05,0.16,0.25,0.5,0.75,0.84,0.95]
    
; SDSS     
    samples = yanny_readone(idlpath+'sdss_samples.par')
    nsamples = n_elements(samples)
    sdssstats = replicate(stats,nsamples)
    for ii = 0, nsamples-1 do begin
       these = where((sdsskcorr.z ge samples[ii].zmin) and $
         (sdsskcorr.z le samples[ii].zmax) and $
         (sdsskcorr.ugriz_absmag[1] le samples[ii].faintcut) and $
         (sdsskcorr.ugriz_absmag[1] ge samples[ii].brightcut) and $
         (sdssispec.dens_edge eq 0),ngal)
       sdssstats[ii].z = mean([samples[ii].zmin,samples[ii].zmax])
       sdssstats[ii].zerr = stddev([samples[ii].zmin,samples[ii].zmax])
       sdssstats[ii].ne_mean = djs_mean(sdssispec[these].dens)
       ss = weighted_quantile(sdssispec[these].dens,quant=quant)
       sdssstats[ii].ne_p05 = ss[0]
       sdssstats[ii].ne_p16 = ss[1]
       sdssstats[ii].ne_p25 = ss[2]
       sdssstats[ii].ne_p50 = ss[3]
       sdssstats[ii].ne_p75 = ss[4]
       sdssstats[ii].ne_p84 = ss[5]
       sdssstats[ii].ne_p95 = ss[6]
       sdssstats[ii].ngal = ngal
    endfor
    im_mwrfits, sdssstats, oiinepath+'oiine_stats_sdss.fits'

; DEEP2 
    samples = yanny_readone(idlpath+'deep2_samples.par')
    nsamples = n_elements(samples)
    deepstats = replicate(stats,nsamples)
    for ii = 0, nsamples-1 do begin
       these = where((deepkcorr.z ge samples[ii].zmin) and $
         (deepkcorr.z le samples[ii].zmax) and $
         (deepkcorr.ugriz_absmag[1] le samples[ii].faintcut) and $
         (deepkcorr.ugriz_absmag[1] ge samples[ii].brightcut) and $
         (deepispec.dens_edge eq 0),ngal)
       deepstats[ii].z = mean([samples[ii].zmin,samples[ii].zmax])
       deepstats[ii].zerr = stddev([samples[ii].zmin,samples[ii].zmax])
       deepstats[ii].ne_mean = djs_mean(deepispec[these].dens)
       ss = weighted_quantile(deepispec[these].dens,quant=quant)
       deepstats[ii].ne_p05 = ss[0]
       deepstats[ii].ne_p16 = ss[1]
       deepstats[ii].ne_p25 = ss[2]
       deepstats[ii].ne_p50 = ss[3]
       deepstats[ii].ne_p75 = ss[4]
       deepstats[ii].ne_p84 = ss[5]
       deepstats[ii].ne_p95 = ss[6]
       deepstats[ii].ngal = ngal
;      if ii eq 2 then stop
    endfor
    im_mwrfits, deepstats, oiinepath+'oiine_stats_deep.fits'

return
end
    
