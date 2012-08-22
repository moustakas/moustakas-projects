pro build_z11_supergrid, supergrid, make_montegrid=make_montegrid, clobber=clobber
; jm12aug14siena - build all the SFH grids we are going to need using
; the supergrid parameter file

    isedfit_sfhgrid_dir = clash_path(/z11)+'montegrids/'
    sfhgrid_paramfile = getenv('CLASH_DIR')+'/z11/z11_sfhgrid.par'
    supergrid_paramfile = getenv('CLASH_DIR')+'/z11/z11_supergrid.par'

    super = yanny_readone(supergrid_paramfile)
    nsuper = n_elements(super)
    struct_print, super

    for ii = 0, nsuper-1 do begin
       build_isedfit_sfhgrid, super[ii].sfhgrid, synthmodels=strtrim(super[ii].synthmodels,2), $
         imf=strtrim(super[ii].imf,2), redcurve=super[ii].redcurve, make_montegrid=make_montegrid, $
         clobber=clobber, isedfit_sfhgrid_dir=isedfit_sfhgrid_dir, sfhgrid_paramfile=sfhgrid_paramfile
    endfor

return
end
