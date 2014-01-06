function bcgsfhs_path, propath=propath, paper=paper, bcg=bcg, skysub=skysub, $
  sersic=sersic, ellipse=ellipse, colormosaics=colormosaics, isedfit=isedfit, $
  ancillary=ancillary
; jm13oct19siena 

    datapath = getenv('IM_ARCHIVE_DIR')+'/projects/clash/bcgsfhs/'
    if keyword_set(skysub) then path = datapath+'skysub/'
    if keyword_set(bcg) then path = datapath+'bcg/'

    path = getenv('IM_PROJECTS_DIR')+'/clash/bcgsfhs/' ; default
    if keyword_set(ellipse) then path = path+'ellipse/'
    if keyword_set(sersic) then path = path+'sersic/'
    if keyword_set(colormosaics) then path = path+'colormosaics/'
    if keyword_set(isedfit) then path = path+'isedfit/'
    if keyword_set(ancillary) then path = path+'ancillary/'

    if keyword_set(propath) then path = getenv('CLASH_DIR')+'/bcgsfhs/'
    if keyword_set(paper) then path = getenv('IM_PAPERS_DIR')+'/projects/clash/bcgsfhs/'
    
return, path
end
