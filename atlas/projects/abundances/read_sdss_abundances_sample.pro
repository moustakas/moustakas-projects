function read_sdss_abundances_sample, sdssnodust=sdssnodust, sdssancillary=sdssancillary, agn=agn, all=all
; jm06apr16uofa - based on READ_SDSS_MZ_SAMPLEn

    datapath = atlas_path(/projects)+'abundances/'
    
    if keyword_set(all) then begin

       file_hii = 'sdss_abundances_hii_speclinefit.fits.gz'
       filenodust_hii = 'sdss_abundances_hii_speclinefit_nodust.fits.gz'
       ancillaryfile_hii = 'sdss_abundances_hii_ancillary.fits.gz'

       file_agn = 'sdss_abundances_agn_speclinefit.fits.gz' 
       filenodust_agn = 'sdss_abundances_agn_speclinefit_nodust.fits.gz'
       ancillaryfile_agn = 'sdss_abundances_agn_ancillary.fits.gz'

       sdssdust_hii = mrdfits(datapath+file_hii,1,/silent)
       sdssdust_agn = mrdfits(datapath+file_agn,1,/silent)
       sdssdust = struct_append(temporary(sdssdust_hii),temporary(sdssdust_agn))

       if arg_present(sdssnodust) then begin
          sdssnodust_hii = mrdfits(datapath+filenodust_hii,1,/silent)
          sdssnodust_agn = mrdfits(datapath+filenodust_agn,1,/silent)
          sdssnodust = struct_append(temporary(sdssnodust_hii),temporary(sdssnodust_agn))
       endif

       if arg_present(sdssancillary) then begin
          sdssancillary_hii = mrdfits(datapath+ancillaryfile_hii,1,/silent)
          sdssancillary_agn = mrdfits(datapath+ancillaryfile_agn,1,/silent)
          sdssancillary = struct_append(temporary(sdssancillary_hii),temporary(sdssancillary_agn))
       endif

       return, sdssdust

    endif

    if keyword_set(agn) then begin
       file = 'sdss_abundances_agn_speclinefit.fits.gz' 
       filenodust = 'sdss_abundances_agn_speclinefit_nodust.fits.gz'
       ancillaryfile = 'sdss_abundances_agn_ancillary.fits.gz'
    endif else begin
       file = 'sdss_abundances_hii_speclinefit.fits.gz'
       filenodust = 'sdss_abundances_hii_speclinefit_nodust.fits.gz'
       ancillaryfile = 'sdss_abundances_hii_ancillary.fits.gz'
    endelse

    sdssdust = mrdfits(datapath+file,1,/silent)
    if arg_present(sdssnodust) then sdssnodust = mrdfits(datapath+filenodust,1,/silent)
    if arg_present(sdssancillary) then sdssancillary = mrdfits(datapath+ancillaryfile,1,/silent)

return, sdssdust
end
    
