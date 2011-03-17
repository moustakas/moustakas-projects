pro readfits_control
; jm01jan13uofa
; read in and display a fits file for the user with ATV

    fitsid = fsc_fileselect(sirtf.widget_ids.base_id,filter='*.fits',/read,$
                            /mustexist,/nomaxsize,objectref=object)

; this won't work for now
    
    object->getfilename=fitsfile
    atv, fitsfile    

return
end
