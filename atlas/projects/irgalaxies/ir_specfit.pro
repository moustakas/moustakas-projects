pro ir_specfit
; jm05aug15uofa - fit the IR galaxies with 4 templates

    datapath = atlas_path(/projects)+'irgalaxies/spec1d/'
    specfitpath = atlas_path(/projects)+'irgalaxies/specfit/'
    ir = read_irgalaxies()
    
; the following code generated the templates used 

;   write_bc03_templates, agegrid=[10.0,100.0,1000.0,10000.0], $
;     tpath=specfitpath, suffix='04', /write, /salpeter

    eigendir = specfitpath
    eigenfile = 'BC03_Z02_salpeter_04_templates.fits'

; fit the IR galaxies    

    suffix = 'ir_04'
    speclist = strtrim(ir.drift_file,2)

    specdata = ispeclinefit(speclist,specres=8.0,snrcut=1.0,dustmodel=0,$
      datapath=datapath,linepath=specfitpath,suffix=suffix,Zmulti=Zmulti,$
      /charlot,/zcrosscor,/postscript,/write,vmaxshift=1000.0,/nologfile,$
      starvdisp=100.0,eigendir=eigendir,eigenfile=eigenfile)

    ir_parse_specfit, datapath=specfitpath
       
stop    
    
return
end    
