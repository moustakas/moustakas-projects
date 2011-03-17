pro write_sings_bc03_templates, bcinfo, ngc1705=ngc1705
; jm06feb12uofa - write out population synthesis templates for AGES

    version = sings_version(/templates)
    
    if keyword_set(ngc1705) then begin

       agegrid = [5.0,10.0,12.0,15.0,20.0]
       suffix = 'ngc1705_'+version

       write_bc03_templates, bc, bcwave, bcinfo, agegrid=agegrid, $
         suffix=suffix, /salpeter, tpath=sings_path(/specfit), $
         minwave=3000.0, maxwave=9000.0, /write

    endif else begin

;      agegrid = [5.0,25.0,100.0,250.0,640.0,1500.0,12000.0]
       agegrid = [5.0,25.0,50.0,125.0,300.0,650.0,1500.0,3000.0,7000.0,13000.0]
       suffix = 'sings_'+version

       write_bc03_templates, bc, bcwave, bcinfo, agegrid=agegrid, $
         suffix=suffix, /salpeter, tpath=sings_path(/specfit), $
         minwave=3000.0, maxwave=9000.0, /write

    endelse
    
return
end
    
