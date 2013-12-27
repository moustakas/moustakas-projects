pro deep2_merge_gandalf_specfit_dr4, specdata, clobber=clobber
; jm13dec20siena - written

    version = deep2_version(/ppxf)
    specfitpath = deep2_path(/ppxf,/dr4)
    spec1dpath = deep2_path(/dr4)
    outpath = deep2_path(/catalogs)

    specdata_outfile = outpath+'deep2.ppxf.specdata.'+version+'.fits'
    specdatafiles = file_search(specfitpath+'specdata_????.fits.gz',count=fcount)

    for ii = 0, fcount-1 do begin
       print, format='("Merging mask ",I0,"/",I0,A10,$)', ii+1, $
         fcount, string(13b)
       specdata1 = mrdfits(specdatafiles[ii],1,/silent)
       if (ii eq 0) then specdata = specdata1 else $
         specdata = [temporary(specdata),specdata1]
    endfor

    im_mwrfits, specdata, specdata_outfile, /clobber
;   im_mwrfits, specfit, specfit_outfile, /create

return    
end
