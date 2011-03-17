function nfgs_path, analysis=analysis, spec1d=spec1d, specfit=specfit, $
  original_spec1d=original_spec1d, ned=ned
; jm05jul24uofa

    home = getenv('RESEARCHPATH')+'/projects/atlas'
    path = home+'/projects/nfgs/'

;   home = getenv('HOME')
;   path = home+'/kennicutt/projects/nfgs/'

    if keyword_set(spec1d) then path = path+'spec1d/'
    if keyword_set(analysis) then path = path+'analysis/'
    if keyword_set(ned) then path = path+'analysis/ned/'
    if keyword_set(specfit) then path = path+'specfit/'
    if keyword_set(original_spec1d) then path = path+'original_spec1d/'
    
return, path
end
