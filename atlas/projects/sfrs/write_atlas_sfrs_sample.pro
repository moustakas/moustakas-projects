pro write_atlas_sfrs_sample, atlasdust, atlasnodust, write=write
; jm05aug26uofa - write an ATLAS sample for the SFRS paper

    outpath = atlas_path(/projects)+'sfrs/'

    sigmacut = 3.0
    
    if (n_elements(atlasdust) eq 0L) then atlasdust = read_integrated(atlasnodust=atlasnodust)

    hakeep = where(atlasdust.h_alpha[0]/atlasdust.h_alpha[1] gt sigmacut,nha)
    hbkeep = where(atlasdust.h_beta[0]/atlasdust.h_beta[1] gt sigmacut,nhb)

    keep = cmset_op(hakeep,'and',hbkeep)
    ngalaxy = n_elements(keep)

    splog, string(ngalaxy,format='(I0)')+'/'+string(n_elements(atlasdust),format='(I0)')+$
      ' galaxies made the S/N cut.'

    hii = where(strtrim(atlasdust[keep].bpt_class,2) eq 'HII',nhii)
    agn = where(strtrim(atlasdust[keep].bpt_class,2) eq 'AGN',nagn)

    splog, 'Found '+string(nhii,format='(I0)')+'/'+string(ngalaxy,format='(I0)')+$
      ' star-forming galaxies.'
    splog, 'Found '+string(nagn,format='(I0)')+'/'+string(ngalaxy,format='(I0)')+$
      ' AGN.'
    
    if keyword_set(write) then begin
       splog, 'Writing '+outpath+'atlas_sfrs_hii_speclinefit.fits.gz'
       mwrfits, atlasdust[keep[hii]], outpath+'atlas_sfrs_hii_speclinefit.fits', /create
       spawn, ['gzip -f '+outpath+'atlas_sfrs_hii_speclinefit.fits'], /sh

       splog, 'Writing '+outpath+'atlas_sfrs_hii_speclinefit_nodust.fits.gz'
       mwrfits, atlasnodust[keep[hii]], outpath+'atlas_sfrs_hii_speclinefit_nodust.fits', /create
       spawn, ['gzip -f '+outpath+'atlas_sfrs_hii_speclinefit_nodust.fits'], /sh

       splog, 'Writing '+outpath+'atlas_sfrs_agn_speclinefit.fits.gz'
       mwrfits, atlasdust[keep[agn]], outpath+'atlas_sfrs_agn_speclinefit.fits', /create
       spawn, ['gzip -f '+outpath+'atlas_sfrs_agn_speclinefit.fits'], /sh

       splog, 'Writing '+outpath+'atlas_sfrs_agn_speclinefit_nodust.fits.gz'
       mwrfits, atlasnodust[keep[agn]], outpath+'atlas_sfrs_agn_speclinefit_nodust.fits', /create
       spawn, ['gzip -f '+outpath+'atlas_sfrs_agn_speclinefit_nodust.fits'], /sh
    endif

return
end
    
