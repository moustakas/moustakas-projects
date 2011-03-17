function gdds_path, spec1d=spec1d, analysis=analysis, raw=raw
; jm05feb08uofa

    home = getenv('HOME')

    rootpath = home+'/research/projects/GDDS/'
    datapath = rootpath
    if keyword_set(spec1d) then datapath = rootpath+'spec1d/'
    if keyword_set(raw) then datapath = rootpath+'raw/'
    if keyword_set(analysis) then datapath = rootpath+'analysis/'
    
return, datapath
end
