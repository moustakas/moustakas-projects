function streams_path, propath=propath, paper=paper, skysub=skysub, $
  psfs=psfs, bcg=bcg, dimage=dimage, postman_bcg=postman_bcg
; jm13jun10siena

    path = getenv('IM_ARCHIVE_DIR')+'/projects/clash/streams/' ; default
    if keyword_set(skysub) then path = path+'skysub/'
    if keyword_set(psfs) then path = path+'psfs/'
    if keyword_set(dimage) then path = path+'dimage/'
    if keyword_set(bcg) then path = path+'bcg/'
    if keyword_set(postman_bcg) then path = path+'bcg_postman/'

    if keyword_set(propath) then path = getenv('CLASH_DIR')+'/streams/'
    if keyword_set(paper) then path = getenv('IM_PAPERS_DIR')+'/projects/clash/streams/'
    
return, path
end
