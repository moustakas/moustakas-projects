pro deep2_merge_gandalf_specfit_dr4, specdata, fixoii=fixoii, clobber=clobber
; jm13dec20siena - written

    version = deep2_version(/ppxf)
    specfitpath = deep2_path(/ppxf,/dr4)
    spec1dpath = deep2_path(/dr4)
    outpath = deep2_path(/catalogs)

    if keyword_set(fixoii) then begin
       specdata_outfile = outpath+'deep2.ppxf.specdata.fixoii.'+version+'.fits'
       specdatafiles = file_search(specfitpath+'specdata_fixoii_????.fits.gz',count=fcount)
    endif else begin
       specdata_outfile = outpath+'deep2.ppxf.specdata.'+version+'.fits'
       specdatafiles = file_search(specfitpath+'specdata_????.fits.gz',count=fcount)
    endelse

    for ii = 0, fcount-1 do begin
       print, format='("Merging mask ",I0,"/",I0,A10,$)', ii+1, $
         fcount, string(13b)
       specdata1 = mrdfits(specdatafiles[ii],1,/silent)
       if (ii eq 0) then specdata = specdata1 else $
         specdata = [temporary(specdata),specdata1]
    endfor

;; if /FIXOII, replace [OII] with the version of the catalog where
;; [OII] is allowed to vary, except where one component of the 
;    if keyword_set(fixoii) then begin
;       specdata1 = mrdfits(outpath+'deep2.ppxf.specdata.'+version+'.fits.gz',1)
;       these = where(specdata1.oii_3727_1_amp[1] gt 0 and specdata1.oii_3727_2_amp[1] gt 0,ngood)
;    endif
    
    im_mwrfits, specdata, specdata_outfile, /clobber
;   im_mwrfits, specfit, specfit_outfile, /create

      
    
    
return    
end
