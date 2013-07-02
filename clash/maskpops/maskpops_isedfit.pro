pro maskpops_isedfit, prelim=prelim, models=models, isedfit=isedfit, $
  qaplot=qaplot, clobber=clobber, thissupergrid=thissupergrid
; jm12may08ucsd
; jm13jul01siena - updated to latest iSEDfit

    prefix = 'maskpops'
    isedfit_dir = maskpops_path(/isedfit)
    montegrids_dir = maskpops_path(/montegrids)
    isedfit_paramfile = isedfit_dir+prefix+'_paramfile.par'
    sfhgrid_paramfile = isedfit_dir+prefix+'_sfhgrid.par'
    supergrid_paramfile = isedfit_dir+prefix+'_supergrid.par'
    
; gather the photometry
    cat = read_maskpops()
;   use_redshift = cat.z ; custom redshift array
    zmin = fix(min(cat.z*10))/10.0
    zmax = ceil(max(cat.z*10))/10.0
    nzz = 3
    zlog = 0

    filterlist = clash_filterlist()
    nfilt = n_elements(filterlist)
    
; supergrid priors
    sfhgrid = 1
    supergrid = 1
    synthmodels = 'fsps'
    imf = 'chab'
    redcurve = 0 ; Calzetti
    igm = 1
    
; SFHgrid priors
    nage = 50 ; 20
    nmonte = 1000 ; 500
    Z = [0.0002,0.03]
    minage = 0.05
    maxage = 6.0
    AV = [0.0,3.0]
    flatAV = 1
    delay_tau = [0.01,6.0]

; --------------------------------------------------
; do the preliminaries: build the parameter files and the Monte Carlo
; grids
    if keyword_set(prelim) then begin
       write_isedfit_paramfile, filterlist, prefix=prefix, minz=zmin, $
         maxz=zmax, nzz=nzz, igm=igm, isedfit_dir=isedfit_dir, clobber=clobber

; delayed
       write_sfhgrid_paramfile, sfhgrid_paramfile, params, clobber=clobber, $
         sfhgrid=1, nage=nage, nmonte=nmonte, Z=Z, minage=minage, $
         maxage=maxage, AV=AV, tau=delay_tau, flatAV=flatAV, /delayed

       write_supergrid_paramfile, supergrid_paramfile, supergrid=supergrid, $
         sfhgrid=sfhgrid, synthmodels=synthmodels, imf=imf, redcurve=redcurve, $
         clobber=clobber

       build_montegrids, sfhgrid_paramfile, supergrid_paramfile=supergrid_paramfile, $
         isedfit_dir=isedfit_dir, clobber=clobber
    endif

; --------------------------------------------------
; build the models
    if keyword_set(models) then begin
       isedfit_models, isedfit_paramfile, isedfit_dir=isedfit_dir, $
         supergrid_paramfile=supergrid_paramfile, clobber=clobber, $
         use_redshift=use_redshift, thissupergrid=thissupergrid
    endif

; --------------------------------------------------
; do the fitting!
    if keyword_set(isedfit) then begin
       maskpops_to_maggies, cat, maggies, ivarmaggies
       isedfit, isedfit_paramfile, maggies, ivarmaggies, cat.z, $
         result, isedfit_dir=isedfit_dir, isedfit_outfile=isedfit_outfile, $
         supergrid_paramfile=supergrid_paramfile, clobber=clobber, $
         sfhgrid_paramfile=sfhgrid_paramfile, outprefix=outprefix, $
         thissupergrid=thissupergrid
    endif 

; --------------------------------------------------
; make a QAplot
    if keyword_set(qaplot) then begin
       yrange = [30,19]
       isedfit_qaplot, isedfit_paramfile, supergrid_paramfile=supergrid_paramfile, $
         isedfit_dir=isedfit_dir, montegrids_dir=montegrids_dir, $
         galaxy=cat.prefix, outprefix=outprefix, clobber=clobber, $
         yrange=yrange, /xlog, thissupergrid=thissupergrid
    endif
    
return
end
