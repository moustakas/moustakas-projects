pro build_maskpops_supergrid, supergrid, make_montegrid=make_montegrid, clobber=clobber
; jm12may08ucsd - build all the SFH grids we are going to need using
; the supergrid parameter file

    isedfit_sfhgrid_dir = maskpops_path(/montegrids)
    sfhgrid_paramfile = getenv('CLASH_DIR')+'/maskpops/maskpops_sfhgrid.par'

    super = get_maskpops_supergrid(supergrid,nsuper=nsuper)
    struct_print, super

    for ii = 0, nsuper-1 do begin
       build_isedfit_sfhgrid, super[ii].sfhgrid, synthmodels=strtrim(super[ii].synthmodels,2), $
         imf=strtrim(super[ii].imf,2), redcurve=super[ii].redcurve, make_montegrid=make_montegrid, $
         clobber=clobber, isedfit_sfhgrid_dir=isedfit_sfhgrid_dir, sfhgrid_paramfile=sfhgrid_paramfile
    endfor

return
end
