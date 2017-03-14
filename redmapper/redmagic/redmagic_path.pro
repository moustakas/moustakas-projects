function redmagic_path, isedfit=isedfit
; jm17feb19siena
    project_root = getenv('IM_ARCHIVE_DIR')+'/projects/redmapper/redmagic/'
    if keyword_set(isedfit) then return, project_root
return, project_root
end
