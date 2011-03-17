function read_nuclear, atlasnodust=atlasnodust, rasort=rasort, silent=silent
; jm03jul6uofa - updated
; jom04mar14uofa - superseded by IRDSPECATLASDUST()
; jm04nov23uofa - impose S/N cuts
; jm05aug26uofa - totally revamped

    version = atlas_version(/specfit)
    path = atlas_path(/spectral_atlas)

    speclinefile = 'nuclear_atlas_speclinefit_'+version+'.fits.gz'
    speclinefilenodust = 'nuclear_atlas_speclinefit_'+version+'_nodust.fits.gz'

    if (not keyword_set(silent)) then splog, 'Reading '+path+speclinefile+'.'
    atlasdust = mrdfits(path+speclinefile,1,/silent)
    if arg_present(atlasnodust) then begin
       if (not keyword_set(silent)) then splog, 'Reading '+path+speclinefilenodust+'.'
       atlasnodust = mrdfits(path+speclinefilenodust,1,/silent)
    endif

return, atlasdust
end
