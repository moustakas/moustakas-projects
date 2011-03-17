function read_calibsfr_sample, silent=silent
; jm10feb21ucsd - read the output from BUILD_CALIBSFR_SAMPLE

    common calibsfr_parent, ages_parent

    datapath = ages_path(/projects)+'calibsfr/'
    thisfile = datapath+'ages_calibsfr.fits.gz'
    if (size(ages_parent,/type) ne 8L) then begin
       if (not keyword_set(silent)) then splog, 'Reading '+thisfile
       ages_parent = mrdfits(thisfile,1,silent=0)
    endif else if (not keyword_set(silent)) then splog, 'Restoring '+file_basename(thisfile)

return, ages_parent
end
