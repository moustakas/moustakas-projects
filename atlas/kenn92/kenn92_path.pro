function kenn92_path, analysis=analysis, ned=ned, spec1d=spec1d, specfit=specfit
; jm05aug02uofa

    home = getenv('HOME')
    
    path = home+'/kennicutt/projects/kenn92'
    if keyword_set(analysis) then path = path+'/analysis/'
    if keyword_set(ned) then path = path+'/analysis/ned/'
    if keyword_set(spec1d) then path = path+'/spec1d/'
    if keyword_set(specfit) then path = path+'/specfit/'
    
return, path
end
