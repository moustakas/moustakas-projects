function bcgsfhs_path, propath=propath, paper=paper, skysub=skysub, $
  bcg=bcg
; jm13oct19siena 

    path = getenv('IM_PROJECTS_DIR')+'/clash/bcgsfhs/' ; default
    datapath = getenv('IM_ARCHIVE_DIR')+'/projects/clash/bcgsfhs/' ; default
    if keyword_set(skysub) then path = datapath+'skysub/'
    if keyword_set(bcg) then path = datapath+'bcg/'

    if keyword_set(propath) then path = getenv('CLASH_DIR')+'/bcgsfhs/'
    if keyword_set(paper) then path = getenv('IM_PAPERS_DIR')+'/projects/clash/bcgsfhs/'
    
return, path
end
