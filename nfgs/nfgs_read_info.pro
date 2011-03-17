function nfgs_read_info
; jm05jul27uofa
    
    path = nfgs_path(/analysis)
    file = 'nfgs_info.fits.gz'

    data = mrdfits(path+file,1,/silent)

return, data
end
