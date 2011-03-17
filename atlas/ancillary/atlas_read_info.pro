function atlas_read_info, silent=silent
; jm05jul27uofa
    
    version = atlas_version(/ancillary)

    path = atlas_path(/analysis)
    file = 'atlas1d_info_'+version+'.fits.gz'

    if (not keyword_set(silent)) then splog, 'Reading '+path+file+'.'
    data = mrdfits(path+file,1,/silent)

return, data
end
