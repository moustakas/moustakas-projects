function read_gto_atlas_sample, datanodust=datanodust
; jm06dec09nyu

    version = gto_log12oh_version()
    datapath = gto_path()

    file = 'gto_atlas_speclinefit_'+version+'.fits.gz' 
    filenodust = 'gto_atlas_speclinefit_nodust_'+version+'.fits.gz' 

    splog, 'Reading '+datapath+file+'.'
    datadust = mrdfits(datapath+file,1,/silent)
    if arg_present(datanodust) then begin
       splog, 'Reading '+datapath+filenodust+'.'
       datanodust = mrdfits(datapath+filenodust,1,/silent)
    endif

return, datadust
end
    
