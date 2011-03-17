; ###########################################################################
; RELEGATED BY SINGS_NED_WEBGET!!!
; ###########################################################################

stop

pro sings_parse_ned, data, photodata
; jm04sep01uofa
; jm05mar22uofa - crossmatch with the RC3    
; jm05jul25uofa - updated

; galaxy_list.txt is the list of objects sent to NED.  the template
; for this batch email is in ${IMPRO_DIR}/ned.  the NED output file is
; called SINGS.NED, which is parsed below into a binary FITS table,
; SINGS_NED.FITS.GZ

    nedpath = sings_path(/analysis)
    
    parse_ned_byname, 'sings.ned', data, nedpath=nedpath, outpath=nedpath, $
      inputnedfile='sings_galaxy_list.txt', outfile='sings_ned.fits'

    parse_ned_photometry, 'sings.photometry.ned', photodata, nedpath=nedpath, outpath=nedpath, $
      inputnedfile='sings_galaxy_list.txt', outfile='sings_ned_photo.fits'

return
end
