function read_sc1120_specfit, galaxy, _extra=extra
; jm05jan28uofa

    if (n_elements(galaxy) eq 0L) then begin
       print, 'specfit = read_sc1120_specfit(galaxy,_extra=extra)'
       return, -1L
    endif
    
    root = 'sc1120'
    specfitpath = sc1120_path(/analysis)

    specfit = irdspecfit(galaxy,specfitpath=specfitpath,root=root,$
      objtagname='GALAXY',_extra=extra)
    
return, specfit
end
