function deep2_path, catalogs=catalogs, dr3=dr3, dr4=dr4, specfit=specfit, $
  window=window, projects=projects, papers=papers, alphadata=alphadata, $
  ppxf=ppxf, auxfiles=auxfiles
; jm06aug28uofa

    datapath = getenv('IM_ARCHIVE_DIR')+'/data/deep2/'
    paperspath = getenv('IM_PAPERS_DIR')+'/projects/deep2/'

    if keyword_set(papers) then return, paperspath
    if keyword_set(catalogs) then begin
       return, getenv('IM_DATA_DIR')+'/deep2/catalogs/'
;      return, datapath+'catalogs/'
    endif
    if keyword_set(window) then return, getenv('IM_DATA_DIR')+'/deep2/window/'
    if keyword_set(auxfiles) then return, getenv('IM_DATA_DIR')+'/deep2/auxfiles/'
    if keyword_set(specfit) then begin
       if keyword_set(dr3) then return, datapath+'specfit_dr3/'
       if keyword_set(dr4) then return, datapath+'specfit_dr4/'
    endif
    if keyword_set(ppxf) then begin
       if keyword_set(dr3) then return, datapath+'ppxf_dr3/'
       if keyword_set(dr4) then return, datapath+'ppxf_dr4/'
    endif
    if keyword_set(projects) then return, datapath+'projects/'
    if keyword_set(alphadata) then return, datapath+'data/alpha/'

    if keyword_set(dr4) then return, datapath+'dr4/'

    if keyword_set(dr3) then begin
       if file_test('/mount/moon1/ioannis/deep2/',/directory) then $
         datapath = '/mount/moon1/ioannis/deep2/' else begin
          message, 'DEEP2/DR3 data not available on external drive.'
       endelse 
       return, datapath+'dr3/'
    endif

return, datapath
end
