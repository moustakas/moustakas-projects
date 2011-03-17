function highzea_path, spec2d=spec2d, spec1d=spec1d, analysis=analysis, mass=mass, $
  papers=papers
; jm02apr15uofa

    home = '/d1/ioannis/'
;   home = getenv('HOME')
    path = home+'/highzea/'

    if keyword_set(spec1d) then path = home+'/highzea/spec1d/'
    if keyword_set(spec2d) then path = home+'/highzea/spec2d/'
    if keyword_set(analysis) then path = home+'/highzea/analysis/'
    if keyword_set(mass) then path = home+'/highzea/mass/'
    if keyword_set(papers) then path = home+'/highzea/papers/'
    
return, path
end
