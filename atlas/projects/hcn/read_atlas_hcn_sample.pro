function read_atlas_hcn_sample, datanodust=datanodust
; jm06dec11nyu

    datapath = hcn_path()
    
    file = 'atlas_hcn_speclinefit.fits.gz' 
    filenodust = 'atlas_hcn_speclinefit_nodust.fits.gz' 

    datadust = mrdfits(datapath+file,1,/silent)
    if arg_present(datanodust) then datanodust = mrdfits(datapath+filenodust,1,/silent)

return, datadust
end
    
