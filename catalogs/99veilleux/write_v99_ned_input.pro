pro write_v99_ned_input
; jm03sep17uofa
; send this list of galaxies to NED for accurate positions

    data = mrdfits('veilleux99_class.fits',1,/silent)
    galaxy = 'IRAS'+data.name

    openw, lun, 'veilleux99_ned_input.txt', /get_lun
    niceprintf, lun, galaxy
    free_lun, lun

return
end
