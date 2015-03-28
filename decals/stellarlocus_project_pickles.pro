function stellarlocus_project_pickles, filterlist, pickles=pickles
; jm09jul24ucsd - project specified bandpasses onto the Pickles+98
;   spectra 
    
    if (n_elements(pickles) eq 0) then pickles = $
      mrdfits('pickles.fits',1,/silent)
;     mrdfits(getenv('IM_PROJECTS_DIR')+'/decals/dr1/qaplots/pickles.fits',1,/silent)
    nband = n_elements(filterlist)
    nstar = n_elements(pickles)
    
    lambda = k_lambda_to_edges(pickles[0].wave)
    nstar = n_elements(pickles)
    maggies = fltarr(nband,nstar)
    for i = 0L, nstar-1L do maggies[*,i] = k_project_filters(lambda,$
      pickles[i].flux,filterlist=filterlist,/silent)
    
return, maggies
end
