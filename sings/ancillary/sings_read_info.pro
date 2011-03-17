function sings_read_info
; jm05jul27uofa
    
    path = sings_path(/analysis)
    version = sings_version(/ancillary)
    file = 'sings_info_'+version+'.fits.gz'
    data = mrdfits(path+file,1,/silent)

return, data
end
