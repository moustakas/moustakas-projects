pro nfgs_parse_ned, data, photo
; jm04nov29uofa - taken from ATLAS2D_NEDSCRIPT
; jm05jul24uofa - updated
    
    nedpath = nfgs_path(/ned)
    outpath = nfgs_path(/analysis)
    
    parse_ned_byname, 'nfgs.ned', data, nedpath=nedpath, outpath=outpath, $
      inputnedfile='nfgs_ned_input.txt', outfile='nfgs_ned.fits'

    parse_ned_photometry, 'nfgs.photometry.ned', photo, nedpath=nedpath, outpath=outpath, $
      inputnedfile='nfgs_ned_input.txt', outfile='nfgs_ned_photo.fits'

return
end
