pro build_mz_sfhgrid, sfhgrid, make_montegrid=make_montegrid, $
  imf=imf, synthmodels=synthmodels, redcurve=redcurve, clobber=clobber
; jm10jan28ucsd - build all the SFH grids we are going to need

    isedfit_sfhgrid_dir = mz_path(/monte)
    sfhgrid_paramfile = getenv('IDL_PROJECTS_DIR')+'/ages/projects/mz/mz_sfhgrid.par'

; defaults    
    if (n_elements(imf) eq 0) then imf = 'chab'
    if (n_elements(synthmodels) eq 0) then synthmodels = 'bc03'
    if (n_elements(redcurve) eq 0) then redcurve = 1 ; charlot

    mzgrid = read_sfhgrid_paramfile(sfhgrid,sfhgrid_paramfile=sfhgrid_paramfile)
    ngrid = n_elements(mzgrid)

    for ii = 0, ngrid-1 do begin
       build_isedfit_sfhgrid, mzgrid[ii].sfhgrid, synthmodels=synthmodels, imf=imf, $
         redcurve=redcurve, make_montegrid=make_montegrid, clobber=clobber, $
         isedfit_sfhgrid_dir=isedfit_sfhgrid_dir, sfhgrid_paramfile=sfhgrid_paramfile
    endfor
       
return
end
    
