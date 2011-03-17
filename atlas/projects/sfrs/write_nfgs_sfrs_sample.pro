pro write_nfgs_sfrs_sample, nfgsdust, nfgsnodust, write=write
; jm05aug26uofa - write an NFGS sample for the SFRS paper

    outpath = atlas_path(/projects)+'sfrs/'

    sigmacut = 3.0
    
    if (n_elements(nfgsdust) eq 0L) then nfgsdust = read_nfgs(nfgsnodust=nfgsnodust)

    hakeep = where(nfgsdust.h_alpha[0]/nfgsdust.h_alpha[1] gt sigmacut,nha)
    hbkeep = where(nfgsdust.h_beta[0]/nfgsdust.h_beta[1] gt sigmacut,nhb)

    keep = cmset_op(hakeep,'and',hbkeep)
    ngalaxy = n_elements(keep)

    splog, string(ngalaxy,format='(I0)')+'/'+string(n_elements(nfgsdust),format='(I0)')+$
      ' galaxies made the S/N cut.'

    hii = where(strtrim(nfgsdust[keep].bpt_class,2) eq 'HII',nhii)
    agn = where(strtrim(nfgsdust[keep].bpt_class,2) eq 'AGN',nagn)

    splog, 'Found '+string(nhii,format='(I0)')+'/'+string(ngalaxy,format='(I0)')+$
      ' star-forming galaxies.'
    splog, 'Found '+string(nagn,format='(I0)')+'/'+string(ngalaxy,format='(I0)')+$
      ' AGN.'
    
    if keyword_set(write) then begin
       splog, 'Writing '+outpath+'nfgs_sfrs_hii_speclinefit.fits.gz'
       mwrfits, nfgsdust[keep[hii]], outpath+'nfgs_sfrs_hii_speclinefit.fits', /create
       spawn, ['gzip -f '+outpath+'nfgs_sfrs_hii_speclinefit.fits'], /sh

       splog, 'Writing '+outpath+'nfgs_sfrs_hii_speclinefit_nodust.fits.gz'
       mwrfits, nfgsnodust[keep[hii]], outpath+'nfgs_sfrs_hii_speclinefit_nodust.fits', /create
       spawn, ['gzip -f '+outpath+'nfgs_sfrs_hii_speclinefit_nodust.fits'], /sh

       splog, 'Writing '+outpath+'nfgs_sfrs_agn_speclinefit.fits.gz'
       mwrfits, nfgsdust[keep[agn]], outpath+'nfgs_sfrs_agn_speclinefit.fits', /create
       spawn, ['gzip -f '+outpath+'nfgs_sfrs_agn_speclinefit.fits'], /sh

       splog, 'Writing '+outpath+'nfgs_sfrs_agn_speclinefit_nodust.fits.gz'
       mwrfits, nfgsnodust[keep[agn]], outpath+'nfgs_sfrs_agn_speclinefit_nodust.fits', /create
       spawn, ['gzip -f '+outpath+'nfgs_sfrs_agn_speclinefit_nodust.fits'], /sh
    endif

return
end
    
