function ages_path, papers=papers, specfit=specfit, spec1d=spec1d, analysis=analysis, $
  projects=projects, catalogs=catalogs, mycatalogs=mycatalogs, $
  spectrophotometry=spectrophotometry, rawdata=rawdata, isedfit=isedfit, $
  ppxf=ppxf, thumbs=thumbs, window=window, qaplots=qaplots, psfs=psfs
; jm04jun04uofa

    agespath = getenv('IM_PROJECTS_DIR')+'/ages/'
    paperspath = getenv('IM_PAPERS_DIR')+'/projects/ages/'
    archivepath = getenv('IM_ARCHIVE_DIR')+'/projects/ages/'

; general path names    
    
    if keyword_set(papers) then return, paperspath
    if keyword_set(analysis) then return, agespath+'analysis/'
    if keyword_set(thumbs) then return, agespath+'thumbs/'
    if keyword_set(projects) then return, agespath+'projects/'
    if keyword_set(isedfit) then return, agespath+'isedfit/'
    if keyword_set(catalogs) then return, agespath+'catalogs/'
    if keyword_set(mycatalogs) then return, agespath+'mycatalogs/'
    if keyword_set(window) then return, agespath+'window/'
    if keyword_set(qaplots) then return, agespath+'qaplots/'
    if keyword_set(psfs) then return, agespath+'psfs/'
    if keyword_set(spectcrophotometry) then return, agespath+'spectrophotometry/'

    if keyword_set(ppxf) then return, agespath+'ppxf/'

    if keyword_set(rawdata) then begin
       rawdatapath = archivepath+'rawdata/'
       if (file_test(rawdatapath,/directory) eq 0) then message, $
         'RAWDATA path '+rawdatapath+' does not exist'
       return, rawdatapath
    endif

    if keyword_set(specfit) then begin
       specfitpath = archivepath+'specfit/'
       if (file_test(specfitpath,/directory) eq 0) then message, $
         'SPECFIT path '+specfitpath+' does not exist'
       return, specfitpath
    endif

    if keyword_set(spec1d) then begin
       spec1dpath = archivepath+'spec1d/'
       if (file_test(spec1dpath,/directory) eq 0) then message, $
         'SPEC1D path '+spec1dpath+' does not exist'
       return, spec1dpath
    endif

return, './'    
end
