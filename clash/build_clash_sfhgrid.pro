pro build_clash_sfhgrid, supergrid=supergrid, make_montegrid=make_montegrid, clobber=clobber
; jm11may05ucsd - build all the SFH grids we are going to need using
; the supergrid parameter file

    sfhgrid_basedir = clash_path(/monte)
    sfhgrid_paramfile = clash_path(/ised)+'clash_sfhgrid.par'

    supergrid_paramfile = clash_path(/ised)+'clash_supergrid.par'
    super = yanny_readone(supergrid_paramfile)
    if (n_elements(supergrid) ne 0) then begin
       match2, super.supergrid, supergrid, m1, m2
       if (total(m2 eq -1) ne 0) then message, 'Unknown supergrid!'
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
