function read_04maier, ohcuts=ohcuts
; jm05jan01uofa
    
    root = '04maier'
    
    path = getenv('CATALOGS_DIR')+'/'+root+'/'
    file = root+'.fits.gz'

    data = mrdfits(path+file,1,/silent)
    
    if keyword_set(ohcuts) then begin
       indx = where((data.z_obj gt -900.0) and (data.M_B gt -900.0) and $
         (data.zstrong_r23 gt -900.0) and (data.zstrong_o32 gt -900.0))
       data = data[indx]
    endif
    
return, data
end
