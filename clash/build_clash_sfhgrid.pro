pro build_hizea_sfhgrid, sfhgrid=sfhgrid, make_montegrid=make_montegrid, $
  imf=imf, synthmodels=synthmodels, redcurve=redcurve, clobber=clobber
; jm10jan28ucsd - build all the SFH grids we are going to need

    sfhgrid_basedir = hizea_path(/monte)
    sfhgrid_paramfile = hizea_path(/mass)+'hizea_sfhgrid.par'

; defaults    
    if (n_elements(imf) eq 0) then imf = 'chab'
    if (n_elements(synthmodels) eq 0) then synthmodels = 'bc03'
    if (n_elements(redcurve) eq 0) then redcurve = 1 ; charlot

    hizeagrid = read_sfhgrid_paramfile(sfhgrid,sfhgrid_paramfile=sfhgrid_paramfile)
    ngrid = n_elements(hizeagrid)

    for ii = 0, ngrid-1 do begin
       build_isedfit_sfhgrid, hizeagrid[ii].sfhgrid, synthmodels=synthmodels, imf=imf, $
         redcurve=redcurve, make_montegrid=make_montegrid, clobber=clobber, $
         sfhgrid_basedir=sfhgrid_basedir, sfhgrid_paramfile=sfhgrid_paramfile
    endfor
       
return
end
    
