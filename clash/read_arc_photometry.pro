function read_arc_photometry
; jm11oct14ucsd - read the SE photometry for the lensed galaxies

    catpath = clash_path(/cat)
    sample = rsex(catpath+'arc_sample.dat')
    narc = n_elements(sample)
    
; build a template catalog
    filt = clash_filterlist(short=short)
    nfilt = n_elements(filt)
    cat = {cluster: '', id: 0L, galaxy: '', ra: 0D, dec: 0D, z: 0.0}
    for ii = 0, nfilt-1 do cat = create_struct(cat,$
      short[ii]+'_mag',-99.0,short[ii]+'_magerr',-99.0)
    cat = replicate(cat,narc)
    cat.cluster = sample.cluster

; loop through the clusters in our sample    
    allcluster = strtrim(sample.cluster,2)
    cluster = allcluster[uniq(allcluster,sort(allcluster))]
    ncluster = n_elements(cluster)

;   for ic = ncluster-1, ncluster-1 do begin
    for ic = 0, ncluster-1 do begin
       splog, 'Cluster '+cluster[ic]
       phot = read_clash_catalog(cluster[ic])
       zspec = read_clash_catalog(cluster[ic],/redshift)

       these = where(cluster[ic] eq allcluster,nthese)
       spherematch, phot.ra, phot.dec, sample[these].ra, sample[these].dec, 1D/3600.0, p1, p2
       spherematch, zspec.ra, zspec.dec, sample[these].ra, sample[these].dec, 1D/3600.0, z1, z2

       if (n_elements(p1) ne nthese) or (n_elements(z1) ne nthese) then $
         message, 'Missing photometry and redshifts!'
;      djs_plot, sample[these[p2]].ra, sample[these[p2]].dec, psym=6, sym=3
;      djs_oplot, sample[these[z2]].ra, sample[these[z2]].dec, psym=7, sym=3, color='green'

       cat[these] = im_struct_assign(phot[p1],cat[these],/nozero)
       cat[these] = im_struct_assign(struct_trimtags(zspec[z1],except=['ra','dec','id']),cat[these],/nozero)
    endfor
    cat.galaxy = strtrim(cat.cluster,2)+'-'+string(cat.id,format='(I5.5)')
    
return, cat    
end
