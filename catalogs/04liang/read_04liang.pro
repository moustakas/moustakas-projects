function read_04liang, ohcuts=ohcuts;, datanodust=datanodust
; jm04dec01uofa
    
    root = '04liang'
    
    path = getenv('CATALOGS_DIR')+'/'+root+'/'
    file = root+'.fits.gz'
;   filenodust = root+'_nodust.fits.gz'

    data = mrdfits(path+file,1,/silent)
;   if arg_present(datanodust) then datanodust = mrdfits(path+filenodust,1,/silent)
    
    if keyword_set(ohcuts) then begin
       indx = where((data.z_obj gt -900.0) and (data.m_b gt -900.0) and (data.zstrong_ew_12oh_kk04 gt -900.0))
       data = data[indx]
    endif
    
return, data
end
