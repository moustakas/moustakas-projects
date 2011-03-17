function read_01pettini, datanodust=datanodust
; jm04dec02uofa
    
    root = '01pettini'
    
    path = getenv('CATALOGS_DIR')+'/'+root+'/'
    file = root+'.fits.gz'
    filenodust = root+'_nodust.fits.gz'

    data = mrdfits(path+file,1,/silent)
    if arg_present(datanodust) then datanodust = mrdfits(path+filenodust,1,/silent)
    
return, data
end
