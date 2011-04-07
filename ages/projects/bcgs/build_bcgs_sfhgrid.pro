pro build_bcgs_sfhgrid, supergrid=supergrid, make_montegrid=make_montegrid, clobber=clobber
; jm10jan28ucsd - build all the SFH grids we are going to need using
; the supergrid parameter file (see BCGS_SUPERGRID.PAR)

    sfhgrid_basedir = ages_path(/projects)+'bcgs/isedfit/montegrids/'
    sfhgrid_paramfile = ages_path(/projects)+'bcgs/isedfit/bcgs_sfhgrid.par'

    supergrid_paramfile = ages_path(/projects)+'bcgs/isedfit/bcgs_supergrid.par'
    super = yanny_readone(supergrid_paramfile)
    if (n_elements(supergrid) ne 0) then begin
       match2, super.supergrid, supergrid, m1, m2
       if (total(m2 eq -1) ne 0) then message, 'unknown supergrid!'
       match, super.supergrid, supergrid, m1, m2
       srt = sort(m2) & m1 = m1[srt] & m2 = m2[srt]
       super = super[m1]
    endif
    struct_print, super
    nsuper = n_elements(super)

    for ii = 0, nsuper-1 do begin
       build_isedfit_sfhgrid, super[ii].sfhgrid, synthmodels=strtrim(super[ii].synthmodels,2), $
         imf=strtrim(super[ii].imf,2), redcurve=super[ii].redcurve, make_montegrid=make_montegrid, $
         clobber=clobber, sfhgrid_basedir=sfhgrid_basedir, sfhgrid_paramfile=sfhgrid_paramfile
    endfor

return
end
