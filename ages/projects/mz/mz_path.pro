function mz_path, isedfit=isedfit, montegrids=montegrids, paper=paper, $
  matchedages=matchedages
; jm11apr11ucsd - return various paths

    mzpath = getenv('IM_RESEARCH_DIR')+'/projects/ages/projects/mz/'
    if keyword_set(isedfit) then return, getenv('IM_ARCHIVE_DIR')+$
      '/projects/ages/projects/mz/isedfit/'
    if keyword_set(montegrids) then return, getenv('IM_ARCHIVE_DIR')+$
      '/projects/ages/projects/mz/isedfit/montegrids/'
    if keyword_set(paper) then return, getenv('IM_PAPERS_DIR')+'/projects/ages/mz/FIG_MZ/'

return, mzpath
end
