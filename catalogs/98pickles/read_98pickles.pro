function read_98pickles
;jm08jul25nyu 

    path = filepath('',root_dir=getenv('CATALOGS_DIR'),$
      subdirectory='98pickles')
    file = path+'98pickles.fits'
    
return, mrdfits(file,1,/silent)
end
    
