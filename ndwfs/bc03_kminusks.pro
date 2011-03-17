pro kminusks
; jm09may26nyu - compute the range of K-Ks colors of galaxies     

    filt = ['ndwfs_K.par','twomass_Ks.par']
    v2ab = k_vega2ab(filterlist=filt,/kurucz,/silent)

    sfhpath = getenv('ISEDFIT_SFHGRID_DIR')+'/basemodels/bc03/'
    sfhfiles = 'chab_Z0.02_tau_0'+['0','3','5','7','9']+'.0Gyr.fits.gz'
    nfile = n_elements(sfhfiles)
    kks = fltarr(nfile)
    kks_vega = fltarr(nfile)
    thisage = 12.0 ; [Gyr]
    for ff = 0, nfile-1 do begin
       s1 = mrdfits(sfhpath+sfhfiles[ff],1)
       age = s1.age/1E9
       get_element, age, thisage, indx
       kwave = k_lambda_to_edges(s1.wave)
       m1 = k_project_filters(kwave,s1.flux[*,indx],filterlist=filt)
       m1 = reform(m1)
       kks[ff] = -2.5*alog10(m1[0]/m1[1]) ; K-Ks (AB)
       kks_vega[ff] = kks[ff] + (v2ab[0]-v2ab[1])
    endfor
    niceprint, sfhfiles, kks, kks_vega

return
end
