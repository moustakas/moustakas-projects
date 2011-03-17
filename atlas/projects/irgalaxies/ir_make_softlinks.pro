pro ir_make_softlinks
; jm05aug15uofa - construct softlinks of the IR-selected spectra

    spec1dpath = atlas_path(/atlas1d)
    outpath = atlas_path(/projects)+'irgalaxies/spec1d/'
    
    ir = read_irgalaxies()
    ngalaxy = n_elements(ir)

    for i = 0L, ngalaxy-1L do spawn, ['ln -s '+spec1dpath+strtrim(ir[i].drift_file,2)+' '+$
      outpath+strtrim(ir[i].drift_file,2)], /sh
    
return
end
