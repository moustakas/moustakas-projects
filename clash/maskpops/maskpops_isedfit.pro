pro maskpops_isedfit, write_paramfile=write_paramfile, build_grids=build_grids, $
  models=models, isedfit=isedfit, qaplot=qaplot, clobber=clobber, thissfhgrid=thissfhgrid
; jm12may08ucsd
; jm13jul01siena - updated to latest iSEDfit

    prefix = 'maskpops'
    isedfit_dir = maskpops_path(/isedfit)
    montegrids_dir = maskpops_path(/montegrids)
    isedfit_paramfile = isedfit_dir+prefix+'_paramfile.par'
    
; gather the photometry
    cat = read_maskpops()
    cat.z = cat.z+[-0.01,0.02,-0.03,+0.03,-0.04,0.015,-0.015]
    cat = cat[sort(cat.z)]
    use_redshift = cat.z ; custom redshift array
    
;   zmin = min(use_redshift)
;   zmax = max(use_redshift)
;   nzz = n_elements(use_redshift)
;   zmin = fix(min(cat.z*10))/10.0
;   zmax = ceil(max(cat.z*10))/10.0
;   nzz = 3
;   zlog = 0

    filterlist = maskpops_filterlist()
;   filterlist = clash_filterlist()
    nfilt = n_elements(filterlist)
    
; --------------------------------------------------
; do the preliminaries: build the parameter files and the Monte Carlo
; grids
    if keyword_set(write_paramfile) then begin
       synthmodels = 'pegase'
       imf = 'salp'
;      synthmodels = 'fsps_v2.4_miles'
;      imf = 'chab'
       redcurve = 'calzetti'
       nmodel = 500             ; 30000L
       write_isedfit_paramfile, params=params, isedfit_dir=isedfit_dir, $
         prefix=prefix, filterlist=filterlist, use_redshift=use_redshift, $
         synthmodels=synthmodels, imf=imf, redcurve=redcurve, /igm, $
         sfhgrid=1, nmodel=nmodel, age=[0.05,6.0], tau=[0.01,6.0], $
         Zmetal=[0.0008,0.03], AV=[0.0,3.0], pburst=0.1, interval_pburst=1.0, $
         oiiihb=[-1.0,1.0], /nebular, /flatAV, /delayed, clobber=clobber

       write_isedfit_paramfile, params=params, isedfit_dir=isedfit_dir, $
         prefix=prefix, filterlist=filterlist, use_redshift=use_redshift, $
         synthmodels=synthmodels, imf=imf, redcurve=redcurve, /igm, $
         sfhgrid=2, nmodel=nmodel, age=[0.05,6.0], tau=[0.01,6.0], $
         Zmetal=[0.0008,0.03], AV=[0.0,3.0], pburst=0.1, interval_pburst=1.0, $
         nebular=0, /flatAV, /delayed, /append, clobber=clobber
    endif

    if keyword_set(build_grids) then begin
       isedfit_montegrids, isedfit_paramfile, montegrids_dir=montegrids_dir, $
         thissfhgrid=thissfhgrid, clobber=clobber
    endif

; --------------------------------------------------
; build the models
    if keyword_set(models) then begin
       isedfit_models, isedfit_paramfile, isedfit_dir=isedfit_dir, $
         montegrids_dir=montegrids_dir, clobber=clobber, thissfhgrid=thissfhgrid
    endif

; --------------------------------------------------
; do the fitting!
    if keyword_set(isedfit) then begin
       maskpops_to_maggies, cat, maggies, ivarmaggies
       isedfit, isedfit_paramfile, maggies, ivarmaggies, cat.z, thissfhgrid=thissfhgrid, $
         isedfit_dir=isedfit_dir, outprefix=outprefix, isedfit_results=ised, $
         isedfit_post=isedpost, clobber=clobber
    endif 

; --------------------------------------------------
; make a QAplot
    if keyword_set(qaplot) then begin
       yrange = [30,19]
       isedfit_qaplot, isedfit_paramfile, isedfit_dir=isedfit_dir, $
         montegrids_dir=montegrids_dir, thissfhgrid=thissfhgrid, $
         galaxy=cat.prefix, outprefix=outprefix, clobber=clobber, $
         yrange=yrange, /xlog
    endif
    
return
end
