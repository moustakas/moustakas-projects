function read_ages_gandalf, solar=solar, unfluxed=unfluxed, silent=silent
; jm09nov30ucsd - see AGES_GANDALF_SPECFIT
;   common ages_read_gandalf, ages_ppxf

    version = ages_version(/ppxf_specfit)    
    specfitpath = ages_path(/ppxf)

    if keyword_set(unfluxed) then begin
       thisfile = specfitpath+'ages_specdata_gandalf_unfluxed_'+version+'.fits.gz'
    endif else begin
       if keyword_set(solar) then suffix = '_solar' else suffix = ''
       thisfile = specfitpath+'ages_specdata_gandalf_'+version+suffix+'.fits.gz'
    endelse

    if (size(ppxf,/type) ne 8) then begin
       if (not keyword_set(silent)) then splog, 'Reading '+thisfile
       ppxf = mrdfits(thisfile,1,silent=0)
    endif else if (keyword_set(silent) eq 0) then splog, 'Restoring '+file_basename(thisfile)

return, ppxf    
end
