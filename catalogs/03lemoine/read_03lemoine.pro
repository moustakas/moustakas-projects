function read_03lemoine
; jm05sep10uofa
    
    root = '03lemoine'
    
    path = getenv('CATALOGS_DIR')+'/'+root+'/'
    file = root+'.fits.gz'

    data = mrdfits(path+file,1,/silent)
    
return, data
end
