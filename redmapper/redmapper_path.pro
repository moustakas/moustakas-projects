function redmapper_path, isedfit=isedfit, version=version, $
  paper=paper, redbaryons=redbaryons, icl=icl, qaplots=qaplots, $
  decamlegacy=decamlegacy
; jm13apr08siena
    version = 'v5.10'
    redmapper_path = getenv('REDMAPPER_DATA')+'/'
    paperpath = getenv('IM_PAPERS_DIR')+'/projects/redmapper/'

    if keyword_set(isedfit) then return, getenv('IM_ARCHIVE_DIR')+'/projects/redmapper/'
;   if keyword_set(isedfit) then return, getenv('REDMAPPER_DATA')+'/'

    if keyword_set(paper) then begin
       if keyword_set(icl) then return, paperpath+'icl/'
       if keyword_set(redbaryons) then return, paperpath+'redbaryons/'
    endif

    if keyword_set(qaplots) then begin
       if keyword_set(icl) then return, redmapper_path+'icl/qaplots/'
       if keyword_set(redbaryons) then return, redmapper_path+'redbaryons/qaplots/'
    endif

    if keyword_set(icl) then return, redmapper_path+'icl/'
    if keyword_set(redbaryons) then return, redmapper_path+'redbaryons/'
    if keyword_set(decamlegacy) then return, redmapper_path+'decamlegacy/'

return, redmapper_path
end
