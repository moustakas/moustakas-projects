function read_ediscs_gandalf, solar=solar, silent=silent
; jm10apr28ucsd - see EDISCS_GANDALF_SPECFIT

    version = ediscs_version(/ppxf_specfit)    
    specfitpath = ediscs_path(/ppxf)
    if keyword_set(solar) then suffix = 'solar_' else suffix = ''

    thisfile = specfitpath+'ediscs_specdata_'+suffix+version+'.fits.gz'
    if (size(ppxf,/type) ne 8) then begin
       if (not keyword_set(silent)) then splog, 'Reading '+thisfile
       ppxf = mrdfits(thisfile,1,silent=0)
    endif else if (keyword_set(silent) eq 0) then splog, 'Restoring '+file_basename(thisfile)

return, ppxf
end
