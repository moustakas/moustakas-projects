function sings_path, spec1d=spec1d, spec2d=spec2d, ascii=ascii, analysis=analysis, $
  specfit=specfit, templates=templates, observing=observing, projects=projects, $
  dss=dss, web=web, papers=papers, ppxf=ppxf
; jm03jan16uofa

    datapath = getenv('IM_RESEARCH_DIR')+'/projects/sings/'
    paperspath = getenv('IM_PAPERS_DIR')+'/projects/sings/'

    if keyword_set(spec1d) then return, datapath+'spec1d/' 
    if keyword_set(ascii) then return, datapath+'ascii/' 
    if keyword_set(analysis) then return, datapath+'analysis/' 
    if keyword_set(specfit) then return, datapath+'specfit/' 
    if keyword_set(templates) then return, datapath+'templates/' 
    if keyword_set(observing) then return, datapath+'observing/' 
    if keyword_set(projects) then return, datapath+'projects/' 
    if keyword_set(dss) then return, datapath+'DSS/' 
    if keyword_set(web) then return, datapath+'public_html/sings/' 
    if keyword_set(papers) then return, paperspath
    
    if keyword_set(ppxf) then return, datapath+'ppxf/' 

    if keyword_set(spec2d) then begin
       if file_test('/mount/bias1/ioannis/sings/spec2d/',/directory) then begin
          path = '/mount/bias1/ioannis/sings/spec2d/' 
       endif else begin
          if file_test('/Volumes/WDexternal/data/sings/spec2d/',/directory) then $
            path = '/Volumes/WDexternal/data/sings/spec2d/' else message, 'External drive not mounted!'
       endelse
    endif
    
return, path
end
