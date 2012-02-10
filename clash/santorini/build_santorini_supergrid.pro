pro build_santorini_supergrid, supergrid, make_montegrid=make_montegrid, clobber=clobber
; jm12feb09ucsd - build all the SFH grids we are going to need using
; the supergrid parameter file

    isedfit_sfhgrid_dir = clash_path(/santorini)+'montegrids/'
    sfhgrid_paramfile = getenv('CLASH_DIR')+'/santorini/santorini_sfhgrid.par'

    super = get_santorini_supergrid(supergrid,nsuper=nsuper)
    struct_print, super

    for ii = 0, nsuper-1 do begin
       build_isedfit_sfhgrid, super[ii].sfhgrid, synthmodels=strtrim(super[ii].synthmodels,2), $
         imf=strtrim(super[ii].imf,2), redcurve=super[ii].redcurve, make_montegrid=make_montegrid, $
         clobber=clobber, isedfit_sfhgrid_dir=isedfit_sfhgrid_dir, sfhgrid_paramfile=sfhgrid_paramfile
    endfor

return
end
