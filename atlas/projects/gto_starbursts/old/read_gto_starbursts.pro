function read_gto_starbursts
; jm05jul25uofa
    
    path = gto_path()
    file = 'gto_starbursts_ned.fits.gz'

    data = mrdfits(path+file,1,/silent)

return, data
end
