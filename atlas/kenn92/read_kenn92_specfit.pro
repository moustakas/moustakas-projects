function read_kenn92_specfit, galaxy, _extra=extra
; jm04apr28uofa
    
    specfitpath = atlas_path(/kenn92)+'analysis/'
    root = 'kenn92'

    specfit = irdspecfit(galaxy,specfitpath=specfitpath,root=root,$
      objtagname='GALAXY',_extra=extra)
    
return, specfit
end
