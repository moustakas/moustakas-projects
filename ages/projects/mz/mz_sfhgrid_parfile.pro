function mz_sfhgrid_parfile
; parameter file describing the various grids; used by
; BUILD_MZ_SFHGRID and MZ_ISEDFIT
    parfile = ages_path(/projects)+'mz/mz_sfhgrid.par'
return, parfile
end
