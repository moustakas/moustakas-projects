function ediscs_path, papers=papers, specfit=specfit, $ ; analysis=analysis, $
  projects=projects, spec1d=spec1d, bamford_spec1d=bamford_spec1d, $
  isedfit=isedfit, html=html, sdss=sdss, ppxf=ppxf, catalogs=catalogs, $
  mycatalogs=mycatalogs
; jm06sep27nyu

    datapath = getenv('RESEARCHPATH')+'/projects/ediscs/'
    paperspath = getenv('PAPERSPATH')+'/projects/ediscs/'

    if keyword_set(papers) then return, paperspath
    if keyword_set(specfit) then return, datapath+'specfit/'
;   if keyword_set(analysis) then return, datapath+'analysis/'
    if keyword_set(catalogs) then return, datapath+'catalogs/'
    if keyword_set(mycatalogs) then return, datapath+'mycatalogs/'
    if keyword_set(projects) then return, datapath+'projects/'
    if keyword_set(spec1d) then return, datapath+'spec1d/'
    if keyword_set(bamford_spec1d) then return, datapath+'spec1d/bamford/'
    if keyword_set(isedfit) then return, datapath+'isedfit/'
    if keyword_set(html) then return, datapath+'html/'
    if keyword_set(sdss) then return, datapath+'sdss/'
    
    if keyword_set(ppxf) then return, datapath+'ppxf/'

return, datapath
end
