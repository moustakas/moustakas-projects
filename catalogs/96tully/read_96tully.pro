function read_96tully
; jm05mar01uofa
    
    root = '96tully'
    
    path = getenv('CATALOGS_DIR')+'/'+root+'/'
    file = root+'.fits.gz'

    data = mrdfits(path+file,1,/silent)

return, data
end
