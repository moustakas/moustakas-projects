pro write_ediscs_bc03_templates, bcinfo
; jm07aug30nyu - write out population synthesis templates for EDISCS 
; jm08may20nyu - v1.1

    agegrid = [5.0,25.0,50.0,125.0,300.0,$
      650.0,1500.0,3000.0,7000.0,13000.0] ; v2.0

    suffix = 'ediscs'
    write_bc03_templates, bc, bcwave, bcinfo, agegrid=agegrid, $
      suffix=suffix, tpath=ediscs_path(/specfit), $
      minwave=912.0, maxwave=10000.0, /write

return
end
    
