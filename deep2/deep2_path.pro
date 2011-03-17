function deep2_path, analysis=analysis, dr3=dr3, specfit=specfit, $
  projects=projects, papers=papers, alphadata=alphadata
; jm06aug28uofa

    datapath = getenv('RESEARCHPATH')+'/projects/deep2/'
    paperspath = getenv('PAPERSPATH')+'/projects/deep2/'

    if keyword_set(papers) then return, paperspath
    if keyword_set(analysis) then return, datapath+'analysis/'
    if keyword_set(specfit) then return, datapath+'specfit/'
    if keyword_set(projects) then return, datapath+'projects/'
    if keyword_set(alphadata) then return, datapath+'data/alpha/'

    if keyword_set(dr3) then begin
       if file_test('/mount/moon1/ioannis/deep2/',/directory) then $
         datapath = '/mount/moon1/ioannis/deep2/' else begin
          message, 'DEEP2/DR3 data not available on external drive.'
       endelse 
       return, datapath+'dr3/'
    endif

return, datapath
end
