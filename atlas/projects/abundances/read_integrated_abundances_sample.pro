function read_integrated_abundances_sample, intnodust=intnodust, agn=agn, all=all
; jm05sep05uofa - read the INTEGRATED sample for the sample generated
;                 by WRITE_INTEGRATED_MZ_SAMPLE 

    datapath = atlas_path(/projects)+'abundances/'

    if keyword_set(all) then begin
       splog, 'Not supported right now.'
       return, -1L
;      splog, 'Reading INTEGRATED/All.'
;      file = 'integrated_abundances_speclinefit.fits.gz' 
;      filenodust = 'integrated_abundances_speclinefit_nodust.fits.gz'
    endif

;   if keyword_set(agn) then begin
;      splog, 'Reading INTEGRATED/AGN.'
;      file = 'integrated_abundances_agn_speclinefit.fits.gz' 
;      filenodust = 'integrated_abundances_agn_speclinefit_nodust.fits.gz'
;   endif

    if (not keyword_set(all)) and (not keyword_set(agn)) then begin
       splog, 'Reading INTEGRATED/HII.'
       file = 'integrated_abundances_hii_speclinefit.fits.gz'
       filenodust = 'integrated_abundances_hii_speclinefit_nodust.fits.gz'
    endif

    intdust = mrdfits(datapath+file,1,/silent)
    if arg_present(intnodust) then intnodust = mrdfits(datapath+filenodust,1,/silent)

return, intdust
end
    
