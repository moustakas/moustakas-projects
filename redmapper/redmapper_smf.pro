pro redmapper_smf
; jm13mar29siena - build the stellar mass function for each cluster
    
    prefix = 'redmapper'
    ver = 'v5.2'
    
    catalogs_dir = getenv('REDMAPPER_DATA')+'/catalogs/'
    isedfit_dir = getenv('REDMAPPER_DATA')+'/'
    isedfit_paramfile = isedfit_dir+prefix+'_paramfile.par'
    sfhgrid_paramfile = isedfit_dir+prefix+'_sfhgrid.par'
    supergrid_paramfile = isedfit_dir+prefix+'_supergrid.par'

    thissupergrid = 1
    outprefix = 'testcl1'

    bcgs = mrdfits(catalogs_dir+'dr8_run_redmapper_'+ver+'_lgt20_catalog.fits.gz',1)
    
;   ised = mrdfits(isedfit_dir+'testcl1_fsps_chab_charlot_sfhgrid01.fits.gz',1)
    phot = mrdfits(catalogs_dir+'redmapper_'+ver+'_photometry.fits.gz',1)
    phot = phot[where(phot.mem_match_id lt 20)] ; hack!

    allcl = phot.mem_match_id
    cl = allcl[uniq(allcl,sort(allcl))]
    ncl = n_elements(cl)

    histmin = 8.0
    histmax = 12.5
    binsize = 0.1
    nbin = long((histmax-histmin)/binsize+1)
    massaxis = histmin+lindgen(nbin)*binsize+binsize/2.0

    out = replicate({mem_match_id: 0L, z: 0.0, lambda_chisq: 0.0, $
      lambda_chisq_e: 0.0, nmem: 0, phi: fltarr(nbin)},ncl)
    out.mem_match_id = cl
    
    for ii = 0, 10 do begin
;   for ii = 0, ncl-1 do begin
       index = where(cl[ii] eq allcl,nmem)
       mass = isedfit_reconstruct_posterior(isedfit_paramfile,$
         supergrid_paramfile=supergrid_paramfile,thissupergrid=thissupergrid,$
         isedfit_dir=isedfit_dir,index=index,outprefix=outprefix)
       hist1 = im_hist1d(mass,weight=phot[index].p,histmin=histmin,$
         histmax=histmax,binsize=binsize);,obin=massaxis)
       out[ii].phi = hist1/total(hist1)
       out[ii].nmem = nmem
       out[ii].z = phot[index[0]].z
; match to the BCG catalog to get the richness
       this = where(cl[ii] eq bcgs.mem_match_id)
       out[ii].lambda_chisq = bcgs[this].lambda_chisq
       out[ii].lambda_chisq_e = bcgs[this].lambda_chisq_e
    endfor

    loadct, 3
    djs_plot, massaxis, out[0].phi, psym=10
    for ii = 1, 10 do djs_oplot, massaxis, out[ii].phi, psym=10
    
    
stop    
    
return
end
    
