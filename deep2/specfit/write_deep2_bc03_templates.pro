pro write_deep2_bc03_templates, bcinfo, restrict_old=restrict_old
; jm07sep26nyu - 

    if keyword_set(restrict_old) then begin
       agegrid = [5.0,25.0,100.0,250.0,640.0,1500.0,5000.0]
       suffix = 'restrict_old'
    endif else begin
       agegrid = [5.0,25.0,100.0,250.0,640.0,1500.0,12000.0]
;      suffix = 'deep2'
    endelse

    write_bc03_templates, bc, bcwave, bcinfo, agegrid=agegrid, $
      suffix=suffix, /salpeter, tpath=deep2_path(/specfit), $
      /write

return
end
    
