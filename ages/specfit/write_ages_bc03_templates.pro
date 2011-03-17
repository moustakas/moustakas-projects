pro write_ages_bc03_templates, bcinfo, restrict_old=restrict_old, for_primus=for_primus
; jm06feb12uofa - write out population synthesis templates for AGES

    version = ages_version(/templates)

    if keyword_set(restrict_old) then begin
       agegrid = [5.0,25.0,50.0,125.0,300.0,650.0,1500.0,3000.0,7000.0] ; v2.0
;      agegrid = [5.0,25.0,100.0,250.0,640.0,1500.0,5000.0] ; v1.0
       suffix = 'ages_'+version+'_restrict_old'
    endif else begin
       agegrid = [5.0,25.0,50.0,125.0,300.0,650.0,1500.0,3000.0,7000.0,13000.0] ; v2.0
;      agegrid = [5.0,25.0,100.0,250.0,640.0,1500.0,12000.0] ; v1.0
       suffix = 'ages_'+version
    endelse

    write_bc03_templates, bc, bcwave, bcinfo, agegrid=agegrid, $
      suffix=suffix, /salpeter, tpath=ages_path(/specfit), $
      minwave=912.0, maxwave=25000.0, write=1

return
end
    
