function read_04kobulnicky
; jm05jan01uofa
; jm08apr24nyu - updated
    
    root = '04kobulnicky'
    path = getenv('CATALOGS_DIR')+'/'+root+'/'
    file = root+'.fits'

    data = mrdfits(path+file,1,/silent)
    
return, data
end
