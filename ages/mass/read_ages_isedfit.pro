function read_ages_isedfit, grid=grid, measure=measure, _extra=extra
; jm10feb12ucsd - read the output from AGES_ISEDFIT

    if (n_elements(grid) eq 0) then grid = 2

    isedpath = ages_path(/isedfit)
    paramfile = isedpath+'ages_isedfit.par'
    params = read_isedfit_paramfile(paramfile)
    params.sfhgrid = grid

    fp = isedfit_filepaths(params,iopath=isedpath)
    isedfile = fp.iopath+strtrim(fp.isedfit_outfile,2)+'.gz'
    measurefile = fp.iopath+strtrim(fp.measure_outfile,2)+'.gz'

    splog, 'Reading '+isedfile
    isedfit = mrdfits(isedfile,1,_extra=extra)

    if arg_present(measure) then begin
       splog, 'Reading '+measurefile
       measure = mrdfits(measurefile,1,_extra=extra)
    endif
    
return, isedfit
end
