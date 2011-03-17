function read_89markarian
; jm05jul25uofa
    
    root = '89markarian'
    
    path = getenv('CATALOGS_DIR')+'/'+root+'/'
    file = root+'.fits.gz'

    data = mrdfits(path+file,1,/silent)

return, data
end
