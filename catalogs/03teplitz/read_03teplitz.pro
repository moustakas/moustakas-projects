function read_03teplitz, datanodust=datanodust, samplecuts=samplecuts
; jm04dec01uofa
    
    root = '03teplitz'
    
    path = getenv('CATALOGS_DIR')+'/'+root+'/'
    file = root+'.fits.gz'
    filenodust = root+'_nodust.fits.gz'

    data = mrdfits(path+file,1,/silent)
    if arg_present(datanodust) then datanodust = mrdfits(path+filenodust,1,/silent)

; sample cuts recommended by Teplitz et al. (2003); they also
; recommend requiring S/N > 2.5 in any emission line for a good
; detection 

    if keyword_set(samplecuts) then begin

       good = where(data.z_obj lt 1.415)
       data = data[good]
       if arg_present(datanodust) then datanodust = datanodust[good]

    endif
    
return, data
end
