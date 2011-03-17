pro write_atlas_bc03_templates, bcinfo
; jm07sep27nyu - write out population synthesis templates for AGES

    version = atlas_version(/templates)

;   agegrid = [5.0,25.0,100.0,250.0,640.0,1500.0,12000.0]
    agegrid = [5.0,10.0,25.0,50.0,125.0,300.0,650.0,1500.0,3000.0,7000.0,13000.0]
    suffix = 'atlas_'+version

    write_bc03_templates, bc, bcwave, bcinfo, agegrid=agegrid, $
      suffix=suffix, /salpeter, tpath=atlas_path(/specfit), $
      minwave=912.0, maxwave=7500.0, /write

return
end
    
