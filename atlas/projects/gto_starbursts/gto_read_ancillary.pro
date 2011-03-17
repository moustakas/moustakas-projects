function gto_read_ancillary, silent=silent
; jm06dec12nyu - written
    
    path = gto_path(/ancillary)
    file = 'gto_ancillary_data.fits.gz'

    if (not keyword_set(silent)) then splog, 'Reading '+path+file+'.'
    data = mrdfits(path+file,1,/silent)

return, data
end
