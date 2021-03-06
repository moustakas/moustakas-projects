function read_sdss_sfrs_sample, sdssnodust=sdssnodust, $
  sdssancillary=sdssancillary, agn=agn
; jm05aug24uofa - read the SDSS sample for the SFRs paper generated by
;                 WRITE_SFRS_SDSS_SAMPLE 

    datapath = atlas_path(/projects)+'sfrs/'
    
    if keyword_set(agn) then begin
       file = 'sdss_sfrs_agn_speclinefit.fits.gz' 
       filenodust = 'sdss_sfrs_agn_speclinefit_nodust.fits.gz'
       ancillaryfile = 'sdss_sfrs_agn_ancillary.fits.gz'
    endif else begin
       file = 'sdss_sfrs_hii_speclinefit.fits.gz'
       filenodust = 'sdss_sfrs_hii_speclinefit_nodust.fits.gz'
       ancillaryfile = 'sdss_sfrs_hii_ancillary.fits.gz'
    endelse

    sdssdust = mrdfits(datapath+file,1,/silent)
    if arg_present(sdssnodust) then sdssnodust = mrdfits(datapath+filenodust,1,/silent)
    if arg_present(sdssancillary) then sdssancillary = mrdfits(datapath+ancillaryfile,1,/silent)

return, sdssdust
end
    
