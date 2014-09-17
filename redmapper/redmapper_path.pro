function redmapper_path, isedfit=isedfit, version=version, catalogs=catalogs, $
  paper=paper, redbaryons=redbaryons, icl=icl, decamlegacy=decamlegacy
; jm13apr08siena

    version = 'v5.10'
    project_root = getenv('IM_PROJECTS_DIR')+'/redmapper/'
    archive_root = getenv('IM_ARCHIVE_DIR')+'/projects/redmapper/'
    paperpath = getenv('IM_PAPERS_DIR')+'/projects/redmapper/'
    
    if keyword_set(isedfit) then return, archive_root+version+'/'
    if keyword_set(catalogs) then return, archive_root+'catalogs/'

    if keyword_set(paper) then begin
       if keyword_set(icl) then return, paperpath+'icl/'
       if keyword_set(redbaryons) then return, paperpath+'redbaryons/'
    endif

    if keyword_set(icl) then return, project_root+'icl/'
    if keyword_set(redbaryons) then return, project_root+'redbaryons/'
    if keyword_set(decamlegacy) then return, project_root+'decamlegacy/'

return, project_root
end
