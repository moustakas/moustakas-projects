pro write_v95_ned_input
; jm03sep17uofa
; send this list of galaxies to NED for accurate positions
    
    readcol, 'veilleux95_class.dat', galaxy, class, format='A,A', $
      comment='#', /silent

    openw, lun, 'veilleux95_ned_input.txt', /get_lun
    niceprintf, lun, galaxy
    free_lun, lun

return
end
