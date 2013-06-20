function streams_path, propath=propath, paper=paper
; jm13jun10siena

    path = getenv('IM_ARCHIVE_DIR')+'/projects/clash/streams/' ; default
    if keyword_set(propath) then path = getenv('CLASH_DIR')+'/streams/'
    if keyword_set(paper) then path = getenv('IM_PAPERS_DIR')+'/projects/clash/streams/'
    
return, path
end
