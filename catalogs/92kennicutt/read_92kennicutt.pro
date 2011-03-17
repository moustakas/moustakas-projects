function read_92kennicutt
; jm05mar01uofa

    root = '92kennicutt'
    
    path = getenv('CATALOGS_DIR')+'/'+root+'/'
    file = root+'.fits.gz'

    data = mrdfits(path+file,1,/silent)
    
return, data
end
