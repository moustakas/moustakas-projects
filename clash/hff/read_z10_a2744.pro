function read_z10_a2744, photoz=photoz, bpz_dz=bpz_dz, bpz_redshift=bpz_redshift, $
  index=index
; jm14jun22siena - read the photometry

    path = getenv('CLASH_PROJECTS')+'/hff/z10_a2744/'
    cat = rsex(path+'flx_ape.inp')
    ngal = n_elements(cat)

; add the "best" photometric redshift from iSEDfit 
    zbest = 9.8
    cat = struct_addtags(cat,replicate({galaxy: '', mu: 1.0, z: zbest},ngal))
    cat.galaxy = 'JD1 '+strtrim(cat.id,2)
    cat[3].galaxy = 'JD1A + JD1B + JD1C'

; these magnifications are from Adi
;   cat.mu = [3.87,3.65,1.6,1.0,1.0]
    cat[0:2].mu = [10.01,11.25,3.57]
    
; make a "Total"
;   these = [0,1] ; =A, B
    these = [0,1,2] ; =A, B, C
    tags = tag_names(cat)
    cat1 = im_empty_structure(cat[0],empty_value=-999.0)
;   cat1.galaxy = 'JD1 A+B'
    cat1.galaxy = 'JD1A+JD1B+JD1C'
;   cat1.galaxy = 'JD1 A+B+C'
    cat1.z = zbest
    cat1.mu = 1.0
    for ii = 0, n_elements(tags)-1 do begin
       if strmatch(tags[ii],'*_FLUX') then cat1.(ii) = total(cat[these].(ii))
       if strmatch(tags[ii],'*_FLUXERR') then cat1.(ii) = sqrt(total(cat[these].(ii)^2))
    endfor
;   cat = [cat,cat1]

return, cat
end    
