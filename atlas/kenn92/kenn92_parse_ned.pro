pro kenn92_parse_ned, data, photo
; jm05aug02uofa - taken from NFGS_PARSE_NED
    
    nedpath = kenn92_path(/ned)
    outpath = kenn92_path(/analysis)
    
    parse_ned_byname, 'kenn92.ned', data, nedpath=nedpath, outpath=outpath, $
      inputnedfile='kenn92_ned_input.txt', outfile='kenn92_ned.fits'

    parse_ned_photometry, 'kenn92.photometry.ned', photo, nedpath=nedpath, outpath=outpath, $
      inputnedfile='kenn92_ned_input.txt', outfile='kenn92_ned_photo.fits'

return
end
