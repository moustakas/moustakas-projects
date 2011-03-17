function read_01kobulnicky
; jm05sep10uofa
    
    root = '01kobulnicky'
    
    path = getenv('CATALOGS_DIR')+'/'+root+'/'
    file = root+'.fits.gz'

    data = mrdfits(path+file,1,/silent)
    
return, data
end
