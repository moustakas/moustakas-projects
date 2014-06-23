function read_z10_a2744, photoz=photoz, bpz_dz=bpz_dz, bpz_redshift=bpz_redshift, $
  index=index
; jm14jun22siena - read the photometry

    path = getenv('CLASH_PROJECTS')+'/hff/z10_a2744/'
    cat = rsex(path+'flx_ape.inp')
    ngal = n_elements(cat)

    cat = struct_addtags(cat,replicate({galaxy: '', mu: 1.0, z: 0.0},ngal))
    cat.galaxy = 'JD'+strtrim(cat.id,2)

; I know - you will be needing the magnifications at each point. For
; now you can take 3.87 for image A, 3.65 for image B, and 1.63 for
; image C. That's from Daniel's model.
    cat.mu = [3.87,3.65]
    
; my best redshift    
    cat.z = 9.75
    
return, cat
end    
