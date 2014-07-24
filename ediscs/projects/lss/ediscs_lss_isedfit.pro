pro ediscs_lss_isedfit, write_paramfile=write_paramfile, build_grids=build_grids, $
  model_photometry=model_photometry, isedfit=isedfit, kcorrect=kcorrect, $
  build_oiiflux=build_oiiflux, qaplot_sed=qaplot_sed, thissfhgrid=thissfhgrid, $
  clobber=clobber
; jm14jul22siena - fit the EDisCS sample for Pascale
    
    prefix = 'ediscs_lss'
    isedfit_dir = getenv('IM_PROJECTS_DIR')+'/ediscs/lss/'
    montegrids_dir = isedfit_dir+'montegrids/'
    isedfit_paramfile = isedfit_dir+prefix+'_paramfile.par'

    phot = read_ediscs(/phot)
    spec1d = read_ediscs(/spec1d)
    ngal = n_elements(phot)
    
    ediscs_to_maggies, phot, maggies, ivarmaggies, filterlist=filterlist

    zminmax = [0.05,1.37]
    zbin = 0.05

; --------------------------------------------------
; write the parameter file
    if keyword_set(write_paramfile) then begin
       spsmodels = 'fsps_v2.4_miles'
       imf = 'chab'
       redcurve = 'charlot'
       Zmetal = [0.004,0.03]
       age = [0.1,12.5]
       tau = [0.0,10]
       nmodel = 20000L
       pburst = 0.2
       interval_pburst = 2.0

       write_isedfit_paramfile, params=params, isedfit_dir=isedfit_dir, $
         prefix=prefix, filterlist=filterlist, spsmodels=spsmodels, $
         imf=imf, redcurve=redcurve, igm=0, zminmax=zminmax, nzz=nzz, $
         zbin=zbin, nmodel=nmodel, age=age, tau=tau, Zmetal=Zmetal, $
         pburst=pburst, interval_pburst=interval_pburst, /nebular, /cl
    endif

; --------------------------------------------------
; build the Monte Carlo grids    
    if keyword_set(build_grids) then begin
       isedfit_montegrids, isedfit_paramfile, isedfit_dir=isedfit_dir, $
         montegrids_dir=montegrids_dir, thissfhgrid=thissfhgrid, clobber=clobber;, $
    endif

; --------------------------------------------------
; calculate the model photometry 
    if keyword_set(model_photometry) then begin
       isedfit_models, isedfit_paramfile, isedfit_dir=isedfit_dir, $
         montegrids_dir=montegrids_dir, thissfhgrid=thissfhgrid, $
         clobber=clobber
    endif

; --------------------------------------------------
; fit!
    if keyword_set(isedfit) then begin
       isedfit, isedfit_paramfile, maggies, ivarmaggies, spec1d.z, ra=spec1d.ra, $
         dec=spec1d.dec, isedfit_dir=isedfit_dir, thissfhgrid=thissfhgrid, $
         clobber=clobber, index=index, outprefix=outprefix
    endif 

; --------------------------------------------------
; compute K-corrections
    if keyword_set(kcorrect) then begin
       isedfit_kcorrect, isedfit_paramfile, isedfit_dir=isedfit_dir, $
         montegrids_dir=montegrids_dir, thissfhgrid=thissfhgrid, $
         absmag_filterlist=sdss_filterlist(), band_shift=0.0, $
         clobber=clobber, index=index, outprefix=outprefix
    endif 

; --------------------------------------------------
; generate spectral energy distribution (SED) QAplots
    if keyword_set(qaplot_sed) then begin
       these = shuffle_indx(ngal,num=20)
       isedfit_qaplot_sed, isedfit_paramfile, isedfit_dir=isedfit_dir, $
         montegrids_dir=montegrids_dir, thissfhgrid=thissfhgrid, $
         clobber=clobber, /xlog, galaxy=spec1d.galaxy, index=these
    endif

return
end
