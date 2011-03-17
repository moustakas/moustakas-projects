function read_99veilleux
; jm05may03uofa 

    root = '99veilleux'
    
    path = getenv('CATALOGS_DIR')+'/'+root+'/'
    file = root+'.fits.gz'

    data = mrdfits(path+file,1,/silent)
    
return, data
end
