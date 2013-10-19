function bcgsfhs_path, propath=propath, paper=paper, skysub=skysub
; jm13oct19siena 

    path = getenv('IM_ARCHIVE_DIR')+'/projects/clash/bcgsfhs/' ; default
    if keyword_set(skysub) then path = path+'skysub/'
    if keyword_set(bcg) then path = path+'bcg/'
    if keyword_set(postman_bcg) then path = path+'bcg_postman/'

    if keyword_set(propath) then path = getenv('CLASH_DIR')+'/bcgsfhs/'
    if keyword_set(paper) then path = getenv('IM_PAPERS_DIR')+'/projects/clash/bcgsfhs/'
    
return, path
end
