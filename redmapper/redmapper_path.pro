function redmapper_path, isedfit=isedfit, catalogs=catalogs, version=version
; jm13apr08siena
    redmapper_path = getenv('REDMAPPER_DATA')
    if keyword_set(catalogs) then begin
       version = 'v5.2'
       return, redmapper_path+'/catalogs/'
    endif
    if keyword_set(isedfit) then return, getenv('IM_ARCHIVE_DIR')+'/projects/redmapper/'
;   if keyword_set(isedfit) then return, getenv('REDMAPPER_DATA')+'/'
return, redmapper_path
end
