function read_04pilyugin
; jm05mar11uofa
    
    root = '04pilyugin'
    
    path = getenv('CATALOGS_DIR')+'/'+root+'/'
    file = root+'.fits.gz'

    data = mrdfits(path+file,1,/silent)

return, data
end
