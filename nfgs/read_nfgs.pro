function read_nfgs, nfgsnodust=nfgsnodust, silent=silent
; jm04mar21uofa
; jm04nov23uofa - impose S/N cuts
; jm05aug26uofa - totally revamped
    
    path = getenv('CATALOGS_DIR')+'/nfgs/'
;   path = nfgs_path(/specfit)
    speclinefile = 'nfgs_int_speclinefit.fits.gz'
    speclinefilenodust = 'nfgs_int_speclinefit_nodust.fits.gz'

    if (not keyword_set(silent)) then splog, 'Reading '+path+speclinefile+'.'
    nfgsdust = mrdfits(path+speclinefile,1,/silent)
    if arg_present(nfgsnodust) then begin
       if (not keyword_set(silent)) then splog, 'Reading '+path+speclinefilenodust+'.'
       nfgsnodust = mrdfits(path+speclinefilenodust,1,/silent)
    endif

return, nfgsdust
end
