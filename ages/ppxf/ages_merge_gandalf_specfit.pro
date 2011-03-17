pro ages_merge_gandalf_specfit, specdata, clobber=clobber, $
  solar=solar, unfluxed=unfluxed
; jm09nov30ucsd - written

    version = ages_version(/ppxf_specfit)    
    ppxfpath = ages_path(/ppxf)
    spec1dpath = ages_path(/spec1d)

    if keyword_set(unfluxed) then begin
       specfitpath = spec1dpath+'unfluxed/ppxf/'+version+'/' 
       specdata_outfile = ppxfpath+'ages_specdata_gandalf_unfluxed_'+version+'.fits'
       specdatafiles = file_search(specfitpath+'specdata_???.fits.gz',count=fcount)
    endif else begin
       specfitpath = spec1dpath+'fluxed/tweak/'+version+'/' 
       if keyword_set(solar) then suffix = '_solar' else suffix = ''
       specdata_outfile = ppxfpath+'ages_specdata_gandalf_'+version+suffix+'.fits'
       specdatafiles = file_search(specfitpath+'specdata_???'+suffix+'.fits.gz',count=fcount)
    endelse
    specfit_outfile = repstr(specdata_outfile,'specdata','specfit')
    
    specfitfiles = repstr(specdatafiles,'specdata','specfit')

;   for ii = 0, 5 do begin
    for ii = 0, fcount-1 do begin
       print, format='("Merging plate ",I0,"/",I0,A10,$)', ii+1, $
         fcount, string(13b)
       specdata1 = mrdfits(specdatafiles[ii],1,/silent)
       if (ii eq 0) then specdata = specdata1 else $
         specdata = [temporary(specdata),specdata1]
;      specfit1 = mrdfits(specfitfiles[ii],1,/silent)
;      if (ii eq 0) then specfit = specfit1 else $
;        specfit = [temporary(specfit),specfit1]
    endfor

    im_mwrfits, specdata, specdata_outfile, /clobber
;   im_mwrfits, specfit, specfit_outfile, /create

return    
end
