function read_04davoust
; jm05may20uofa
    
    root = '04davoust'
    
    path = getenv('CATALOGS_DIR')+'/'+root+'/'
    file = root+'_table6.fits.gz'

    data = mrdfits(path+file,1,/silent)

return, data
end
