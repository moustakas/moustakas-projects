pro parse_gto_ned, data, photo
; jm03jul16uofa    
; jm05jul25uofa - updated
; jm06dec07nyu - updated

    nedpath = gto_path(/ancillary)
    outpath = gto_path(/ancillary)

    parse_ned_byname, 'gto.ned', data, nedpath=nedpath, outpath=outpath, $
      inputnedfile='gto_ned_input.txt', outfile='gto_ned.fits'
    parse_ned_photometry, 'gto.photometry.ned', photo, nedpath=nedpath, outpath=outpath, $
      inputnedfile='gto_ned_input.txt', outfile='gto_ned_photo.fits'

return
end
