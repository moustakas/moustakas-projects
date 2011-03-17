pro write_gto_starbursts, atlas, atlas_nodust, write=write
; jm06dec07nyu - write out the subset of ATLAS galaxies that belong to
;                the GTO starburst galaxy subsample

    outpath = gto_path()
    
    if (n_elements(atlas) eq 0L) then atlas = read_integrated(atlasnodust=atlas_nodust)
    indx = index_gto_starbursts(atlas,gto,index_gto=index_gto)

    if keyword_set(write) then begin
       
       
    endif

return
end
