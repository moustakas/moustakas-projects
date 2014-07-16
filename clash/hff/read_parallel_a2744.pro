function read_parallel_a2744, photoz=photoz, bpz_dz=bpz_dz, bpz_redshift=bpz_redshift, $
  index=index
; jm14jul16siena - read the photometry

    date = 'jul16'
    
    path = getenv('CLASH_PROJECTS')+'/hff/parallel_a2744/'
    cat = rsex(path+'a2744p_flx_iso.'+date)
    ngal = n_elements(cat)

; add the "best" photometric redshift from iSEDfit 
    cat = struct_addtags(cat,replicate({galaxy: '', mu: 1.0, z: 0.0},ngal))
    cat.galaxy = 'ID'+strtrim(cat.id,2)

return, cat
end    
