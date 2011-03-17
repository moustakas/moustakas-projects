function read_03gildepaz
; jm05mar01uofa
    
    root = '03gildepaz'
    
    path = getenv('CATALOGS_DIR')+'/'+root+'/'
    file = root+'.fits.gz'

    data = mrdfits(path+file,1,/silent)

return, data
end
