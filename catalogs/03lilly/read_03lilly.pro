function read_03lilly, datanodust=datanodust
; jm05jan01uofa
    
    root = '03lilly'
    
    path = getenv('CATALOGS_DIR')+'/'+root+'/'
    file = root+'.fits'
;   filenodust = root+'_nodust.fits'

    data = mrdfits(path+file,1,/silent)
;   if arg_present(datanodust) then datanodust = mrdfits(path+filenodust,1,/silent)
    
return, data
end
