; ###########################################################################
; RELEGATED BY ATLAS_NED_WEBGET!!!
; ###########################################################################

stop

pro atlas_parse_ned, data, photo
; jm01nov7uofa
; jm02mar22uofa - added some lines 
; jm02may14uofa - added NED photometry 
; jm05jul21uofa - updated

; atlas_ned_input.txt is the list of objects sent to NED.  the
; template for this batch email is in ${IMPRO_DIR}/ned.  the NED
; output file is called ATLAS.NED, which is parsed below into a
; binary FITS table, ATLAS_SAMPLE.FITS

    nedpath = atlas_path(/ned)
    outpath = atlas_path(/analysis)
    
    parse_ned_byname, 'atlas.ned', data, nedpath=nedpath, outpath=outpath, $
      inputnedfile='atlas_ned_input.txt', outfile='atlas_ned.fits'
    parse_ned_photometry, 'atlas.photometry.ned', photo, nedpath=nedpath, outpath=outpath, $
      inputnedfile='atlas_ned_input.txt', outfile='atlas_ned_photo.fits'

; supplement the NED batch results with additional measurements and
; write out again

    info = read_atlas_extra_info()
    data.galaxy = info.galaxy
        
    morez = where((info.redshift ne -1.0) and (data.z lt -900.0),nmorez)
    if (nmorez ne 0L) then begin
       data[morez].z = info[morez].redshift
       splog, 'Additional redshifts for the following objects:'
       niceprint, strtrim(data[morez].galaxy,2)+' ['+strtrim(info[morez].galaxy,2)+']'
    endif

    mored = where((info.dmaj ne -1.0) and (data.dmaj lt -900.0),nmored)
    if (nmored ne 0L) then begin
       data[mored].dmaj = info[mored].dmaj
       data[mored].dmin = info[mored].dmin
       splog, 'Additional diameters for the following objects:'
       niceprint, strtrim(data[mored].galaxy,2)+' ['+strtrim(info[mored].galaxy,2)+']'
    endif

    if (nmorez gt 0L) or (nmored gt 0L) then begin
       mwrfits, data, outpath+'atlas_ned.fits', /create
       spawn, ['gzip -f '+outpath+'atlas_ned.fits'], /sh
    endif

return
end
