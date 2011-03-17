function sg1120_path, papers=papers, analysis=analysis, $
  isedfit=isedfit, mcatalogs=mcatalogs, zcat=zcat, sex=sex
; jm05jan14uofa
; jm07jan17nyu - updated
; jm08jun09nyu - data removed from external drive

    datapath = getenv('RESEARCHPATH')+'/projects/sg1120/'
    paperspath = getenv('PAPERSPATH')+'/projects/sg1120/'

    if keyword_set(papers) then return, paperspath
    if keyword_set(isedfit) then return, datapath+'isedfit/'
    if keyword_set(analysis) then return, datapath+'analysis/'
    if keyword_set(mcatalogs) then return, datapath+'matched_catalogs/'
    if keyword_set(zcat) then return, datapath+'zcat/'
    if keyword_set(sex) then return, datapath+'sex/'
    
return, datapath    
end
