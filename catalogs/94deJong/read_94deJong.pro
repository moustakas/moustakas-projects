function read_94deJong
; jm05mar01uofa
    
    root = '94deJong'
    
    path = getenv('CATALOGS_DIR')+'/'+root+'/'
    file = root+'.fits.gz'

    data = mrdfits(path+file,1,/silent)

return, data
end
