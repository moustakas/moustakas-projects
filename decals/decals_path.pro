function decals_path, dr=dr, tractor=tractor, isedfit=isedfit, $
  gallery=gallery
; jm15mar28siena
    if n_elements(dr) eq 0 then dr = 'dr1'
    
    path = getenv('DECALS_DIR')+'/'+dr+'/'
    if keyword_set(tractor) then path = getenv('DECALS_DIR')+'/'+dr+'/tractor/'
    if keyword_set(gallery) then path = getenv('DECALS_DIR')+'/'+dr+'/gallery-dr1/'
    if keyword_set(isedfit) then path = getenv('DECALS_DIR')+'/isedfit/'+dr+'/'

return, path
end   
