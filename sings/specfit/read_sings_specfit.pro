function read_sings_specfit, galaxy, drift56=drift56, drift20=drift20, $
  nuclear=nuclear, _extra=extra
; jm05jul30uofa
; jm08oct20nyu - updated to the new data model

    version = sings_version(/specfit)
    specfitpath = sings_path(/specfit)+version+'/'

    if (not keyword_set(nuclear)) and (not keyword_set(drift20)) and $
      (not keyword_set(drift56)) then begin
       splog, 'Either NUCLEAR *or* DRIFT20 *or* DRIFT56 keyword must be set!'
       return, -1L
    endif

    if (keyword_set(nuclear) and keyword_set(drift20)) or $
      (keyword_set(nuclear) and keyword_set(drift56)) or $
      (keyword_set(drift20) and keyword_set(drift56)) then begin
       splog, 'Only one keyword (NUCLEAR, DRIFT20, or DRIFT56) can be set at the same time!'
       return, -1L
    endif

    if keyword_set(nuclear) then root = 'sings_nuclear'
    if keyword_set(drift20) then root = 'sings_drift20'
    if keyword_set(drift56) then root = 'sings_drift56'
    if (n_elements(mjdstr) eq 0L) then mjdstr = ''

    allfiles = file_search(djs_filepath(specfitpath+'*'+mjdstr+'_*'+root+'*_specdata.fits.gz'),count=fcount)
    specdatafile = allfiles[(reverse(sort(allfiles)))[0]]
;   splog, 'Reading '+specdatafile
    specdata = mrdfits(specdatafile,1,/silent)

    allfiles = file_search(djs_filepath(specfitpath+'*'+mjdstr+'_*'+root+'*_specfit.fits*'),count=fcount)
    specfitfile = allfiles[(reverse(sort(allfiles)))[0]]
;   splog, 'Reading '+specfitfile
    specfit = mrdfits(specfitfile,0,/silent)

    these = speclinefit_locate(specdata,galaxy)

;   specfit = irdspecfit(galaxy,specfitpath=specfitpath,root=root,$
;     objtagname='GALAXY',_extra=extra)
    
return, reform(specfit[*,*,these])
end
