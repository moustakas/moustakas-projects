function sdss_path, mpa_dr4=mpa_dr4, mpa_dr7=mpa_dr7, $
  projects=projects, isedfit=isedfit, vagc_mpa=vagc_mpa, $
  lowz=lowz
; jm04nov04uofa
; jm05nov02uofa - added DR2, DR4 keywords
; jm09feb23nyu - added MPA_DR7 keyword    

;   datapath = getenv('CATALOGS_DIR')+'/sdss/'
    datapath = getenv('IM_PROJECTS_DIR')+'/sdss/'

    if keyword_set(isedfit) then datapath = getenv('IM_RESEARCH_DIR')+'/projects/sdss/isedfit/'
    if keyword_set(projects) then datapath = getenv('IM_RESEARCH_DIR')+'/projects/sdss/projects/'
    if keyword_set(vagc_mpa) then datapath = getenv('IM_RESEARCH_DIR')+'/data/sdss/vagc-mpa/'
    if keyword_set(lowz) then datapath = getenv('IM_RESEARCH_DIR')+'/data/vagc-lowz-dr6/'

    if keyword_set(mpa_dr7) then datapath = getenv('DATA_DIR')+'/sdss/mpa_dr7_v5_2/'
;   if keyword_set(mpa_dr7) then datapath = getenv('IM_RESEARCH_DIR')+'/data/sdss/mpa_dr7_v5_2/'
;   if keyword_set(mpa_dr4) then datapath = getenv('IM_RESEARCH_DIR')+'/data/sdss/mpa_dr4_v5_1b/'
    
return, datapath
end
