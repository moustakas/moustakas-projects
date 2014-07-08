pro alhambra_isedfit, write_paramfile=write_paramfile, build_grids=build_grids, $
  model_photometry=model_photometry, isedfit=isedfit, kcorrect=kcorrect, $
  build_oiiflux=build_oiiflux, qaplot_sed=qaplot_sed, thissfhgrid=thissfhgrid, $
  clobber=clobber
; jm14jul07siena - fit an ALHAMBRA galaxy sample

; the prefix and ISEDFIT_DIR can be changed arbitrarily    
    prefix = 'alhambra'
    isedfit_dir = getenv('IM_PROJECTS_DIR')+'/alhambra/'
    montegrids_dir = isedfit_dir+'montegrids/'
    isedfit_paramfile = isedfit_dir+prefix+'_paramfile.par'

; read the catalog and convert to maggies
    cat = rsex(isedfit_dir+'sub_SEDs_pop0.dat')
    ngal = n_elements(cat)

    alhambra_to_maggies, cat, maggies, ivarmaggies, filterlist=filterlist

; --------------------------------------------------
; write the parameter file
    if keyword_set(write_paramfile) then begin
       write_isedfit_paramfile, params=params, isedfit_dir=isedfit_dir, $
         prefix=prefix, filterlist=filterlist, spsmodels='fsps_v2.4_miles', $
         imf='chab', redcurve='charlot', igm=0, zminmax=[0.05,1.4], nzz=40, $
         nmodel=1000L, age=[0.1,12.5], tau=[0.0,10], Zmetal=[0.004,0.03], $
         pburst=0.2, interval_pburst=2.0, /nebular, clobber=clobber
    endif

; --------------------------------------------------
; build the Monte Carlo grids    
    if keyword_set(build_grids) then begin
       isedfit_montegrids, isedfit_paramfile, isedfit_dir=isedfit_dir, $
         montegrids_dir=montegrids_dir, thissfhgrid=thissfhgrid, clobber=clobber
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
       isedfit, isedfit_paramfile, maggies, ivarmaggies, $
         cat.z, ra=cat.ra, dec=cat.dec, isedfit_dir=isedfit_dir, $
         thissfhgrid=thissfhgrid, clobber=clobber, index=index, $
         outprefix=outprefix
    endif 

; --------------------------------------------------
; compute K-corrections; specifiy ABSMAG_FILTERLIST() as desired 
    if keyword_set(kcorrect) then begin
       isedfit_kcorrect, isedfit_paramfile, isedfit_dir=isedfit_dir, $
         montegrids_dir=montegrids_dir, thissfhgrid=thissfhgrid, $
         absmag_filterlist=sdss_filterlist(), band_shift=0.0, $
         clobber=clobber, index=index, outprefix=outprefix
    endif 

; --------------------------------------------------
; generate spectral energy distribution (SED) QAplots
    if keyword_set(qaplot_sed) then begin
       galaxy = 'ALH'+strtrim(cat.id,2)
       index = shuffle_indx(ngal,num=25)
       isedfit_qaplot_sed, isedfit_paramfile, isedfit_dir=isedfit_dir, $
         montegrids_dir=montegrids_dir, thissfhgrid=thissfhgrid, $
         clobber=clobber, /xlog, galaxy=galaxy, index=index
    endif

return
end
