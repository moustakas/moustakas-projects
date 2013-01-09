pro z11_mock, supergrid, models=models, isedfit=isedfit, $
  qaplot=qaplot, clobber=clobber
; jm12nov13siena - make a plot

    common com_z11mock, mstar, sfrage, sfr0
    
    
    datapath = clash_path(/z11)
    isedpath = datapath+'isedfit/'
    isedfit_sfhgrid_dir = clash_path(/z11)+'montegrids/'
    paramfile = isedpath+'z11_mock_paramfile.par'
;   sfhgrid_paramfile = getenv('CLASH_DIR')+'/z11/z11_sfhgrid.par'
;   supergrid_paramfile = getenv('CLASH_DIR')+'/z11/z11_supergrid.par'

    supergrid_paramfile = getenv('CLASH_DIR')+'/z11/z11_supergrid.par'
    super = yanny_readone(supergrid_paramfile)

    filters = z11_filterlist(nice=nice_filter)
    ndraw = isedfit_ndraw()

    if n_elements(mstar) eq 0 then begin
       mstar = isedfit_reconstruct_posterior(paramfile,post=post,$
         super=super,isedfit_sfhgrid_dir=isedfit_sfhgrid_dir,isedpath=isedpath,$
         age=age,Z=Z,tau=tau,sfr0=sfr0,av=av,sfrage=sfrage)
    endif

stop    
    

return
end
