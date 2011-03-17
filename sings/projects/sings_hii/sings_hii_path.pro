function sings_hii_path, spec2d=spec2d, spec1d=spec1d, analysis=analysis
; jm06jul31uofa

    home = getenv('HOME')
    path = home+'/kennicutt/sings/projects/sings_hii/'

    if keyword_set(spec1d) then path =   home+'/kennicutt/sings/projects/sings_hii/spec1d/'
    if keyword_set(spec2d) then path =   home+'/kennicutt/sings/projects/sings_hii/spec2d/'
    if keyword_set(analysis) then path = home+'/kennicutt/sings/projects/sings_hii/analysis/'
    
return, path
end
