function read_06mouhcine, ohcuts=ohcuts
; jm05jan01uofa
    
    root = '06mouhcine'
    path = getenv('CATALOGS_DIR')+'/'+root+'/'
    file = root+'.fits.gz'

    data = mrdfits(path+file,1,/silent)
    
    if keyword_set(ohcuts) then begin
       indx = where((data.z_obj gt -900.0) and (data.m_b gt -900.0) and (data.zstrong_ew_12oh_kk04 gt -900.0))
       data = data[indx]
    endif
    
return, data
end
