pro ir_define_sample, atlas, atlasnodust, write=write
; jm05aug12uofa

    if (n_elements(atlas) eq 0L) then atlas = read_integrated(linefitnodust=atlasnodust)

; select the sample    

    ir = where(atlas.ir_lum gt 11.0,nir)
    splog, 'Selecting '+string(nir,format='(I0)')+' IR galaxies.'
    
    iratlas = atlas[ir]
    iratlasnodust = atlasnodust[ir]
    
; write out

    outpath = atlas_path(/projects)+'irgalaxies/'
    outname = 'ir_integrated_speclinefit.fits'
    outnamenodust = 'ir_integrated_speclinefit_nodust.fits'
    
    if keyword_set(write) then begin
    
       splog, 'Writing '+outpath+outname+'.'
       mwrfits, iratlas, outpath+outname, /create
       spawn, ['gzip -f '+outpath+outname], /sh

       splog, 'Writing '+outpath+outnamenodust+'.'
       mwrfits, iratlasnodust, outpath+outnamenodust, /create
       spawn, ['gzip -f '+outpath+outnamenodust], /sh

    endif

return
end
