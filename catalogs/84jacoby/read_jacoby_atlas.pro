function read_jacoby_atlas
; jm03feb3uofa

    jacoby_path = filepath('',root_dir=getenv('CATALOGS_DIR'),$
      subdirectory='84jacoby')

    file = 'jacoby_atlas.fits'
    jacoby = mrdfits(jacoby_path+file,1,/silent)
    
return, jacoby
end    
