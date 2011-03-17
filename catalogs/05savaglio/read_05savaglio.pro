function read_05savaglio, datanodust=datanodust
; jm05sep10uofa
    
    root = '05savaglio'
    
    path = getenv('CATALOGS_DIR')+'/'+root+'/'
    file = root+'.fits'
    filenodust = root+'_nodust.fits'

    data = mrdfits(path+file,1,/silent)
    if arg_present(datanodust) then datanodust = mrdfits(path+filenodust,1,/silent)
    
return, data
end
