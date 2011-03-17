function kenn92_read_info
; jm05aug03uofa
    
    path = kenn92_path(/analysis)
    file = 'kenn92_info.fits.gz'

    data = mrdfits(path+file,1,/silent)

return, data
end
