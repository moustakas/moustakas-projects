function mzr_path, isedfit=isedfit, montegrids=montegrids, paper=paper, $
  matchedages=matchedages
; jm13aug28siena - return various paths

    mzrpath = getenv('IM_RESEARCH_DIR')+'/projects/ages/projects/mzr/'
    if keyword_set(isedfit) then return, getenv('IM_ARCHIVE_DIR')+$
      '/projects/ages/mzr/'
    if keyword_set(montegrids) then return, getenv('IM_ARCHIVE_DIR')+$
      '/projects/ages/mzr/montegrids/'
    if keyword_set(paper) then return, getenv('IM_PAPERS_DIR')+'/projects/ages/mzr/'

return, mzrpath
end
