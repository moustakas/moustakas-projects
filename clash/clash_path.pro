function clash_path, catalogs=catalogs, isedfit=isedfit, montegrids=montegrids
; jm11apr18ucsd - 

    clashpath = getenv('CLASH_DIR')+'/'
    if keyword_set(catalogs) then clashpath = clashpath+'catalogs/'
    if keyword_set(isedfit) then clashpath = clashpath+'isedfit/'
    if keyword_set(montegrids) then clashpath = clashpath+'montegrids/'
    
return, clashpath
end
