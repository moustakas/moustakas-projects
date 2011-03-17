function read_distance_catalog
; jm05feb23uofa
    
    file = filepath('distance_catalog.fits.gz',root=getenv('CATALOGS_DIR'),subdirectory='distances')
    cat = mrdfits(file,1,/silent)
    
return, cat
end
