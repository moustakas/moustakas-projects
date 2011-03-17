function read_05hoyos
; jm06apr20uofa
    
    root = '05hoyos'
    path = getenv('CATALOGS_DIR')+'/'+root+'/'
    file = root+'.fits'

    data = mrdfits(path+file,1,/silent)
    
return, data
end
