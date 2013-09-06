function streams_path, propath=propath, paper=paper, skysub=skysub, $
  psfs=psfs, nobcg=nobcg
; jm13jun10siena

    path = getenv('IM_ARCHIVE_DIR')+'/projects/clash/streams/' ; default
    if keyword_set(skysub) then path = path+'skysub/'
    if keyword_set(psfs) then path = path+'psfs/'
    if keyword_set(nobcg) then path = path+'nobcg/'

    if keyword_set(propath) then path = getenv('CLASH_DIR')+'/streams/'
    if keyword_set(paper) then path = getenv('IM_PAPERS_DIR')+'/projects/clash/streams/'
    
return, path
end
