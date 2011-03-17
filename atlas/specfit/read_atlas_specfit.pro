function read_atlas_specfit, galaxy, nuclear=nuclear, _extra=extra
; jm04mar15uofa
    
    specfitpath = atlas_path(/specfit)

    if keyword_set(nuclear) then $
      root = 'nuclear_atlas' else root = 'integrated_atlas'

    specfit = irdspecfit(galaxy,specfitpath=specfitpath,root=root,$
      objtagname='GALAXY',_extra=extra)
    
return, specfit
end
