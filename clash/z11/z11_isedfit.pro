pro z11_isedfit, write_paramfiles=write_paramfiles, build_montegrids=build_montegrids, $
  models=models, isedfit=isedfit, qaplot=qaplot, clobber=clobber, $
  thissupergrid=thissupergrid
; jm13jul30siena - fit the z=11 candidate in preparation for our
; Spitzer proposal

    prefix = 'test_z11'
    rootpath = getenv('IM_PROJECTS_DIR')+'/clash/z11/'
    isedfit_dir = rootpath
    montegrids_dir = rootpath+'montegrids/'
    isedfit_paramfile = isedfit_dir+prefix+'_paramfile.par'
    sfhgrid_paramfile = isedfit_dir+prefix+'_sfhgrid.par'
    supergrid_paramfile = isedfit_dir+prefix+'_supergrid.par'
    
; gather the photometry
    cat = read_z11()
    zmin = 10.79
    zmax = 10.81
    nzz = 3
    zlog = 0

    filterlist = z11_filterlist()
    nfilt = n_elements(filterlist)
    
; supergrid priors
    sfhgrid = 1
    supergrid = 1
    synthmodels = 'fsps_v2.4_miles'
    igm = 1

    imf = 'chab'
    redcurve = -1 ; 0 ; Calzetti
    
; SFHgrid priors
    nage = 40
    nmonte = 1250
;   Z = [0.004,0.004]
    Z = [0.0002,0.02]
    oiiihb = [-0.2,0.8]
    minage = 0.005
    maxage = 0.45
    AV = [0.0,0.0]
    flatAV = 1
    delay_tau = [0.01,1.0]

; --------------------------------------------------
; build the parameter files
    if keyword_set(write_paramfiles) then begin
       write_isedfit_paramfile, filterlist, prefix=prefix, minz=zmin, $
         maxz=zmax, nzz=nzz, igm=igm, isedfit_dir=isedfit_dir, clobber=clobber

       write_sfhgrid_paramfile, sfhgrid_paramfile, params, clobber=clobber, $
         sfhgrid=1, nage=nage, nmonte=nmonte, Z=Z, minage=minage, $
         maxage=maxage, AV=AV, tau=delay_tau, flatAV=flatAV, /delayed
       write_sfhgrid_paramfile, sfhgrid_paramfile, params, clobber=clobber, $
         sfhgrid=2, nage=nage, nmonte=nmonte, Z=Z, minage=minage, $
         maxage=maxage, AV=AV, tau=delay_tau, flatAV=flatAV, /delayed, $
         oiiihb=oiiihb, /nebular, /append

       write_supergrid_paramfile, supergrid_paramfile, supergrid=[1,2], $
         sfhgrid=[1,2], synthmodels=synthmodels, imf=imf, redcurve=redcurve, $
         clobber=clobber
    endif

; build the Monte Carlo grids
    if keyword_set(build_montegrids) then begin
       build_montegrids, sfhgrid_paramfile, supergrid_paramfile=supergrid_paramfile, $
         isedfit_dir=isedfit_dir, clobber=clobber, thissupergrid=thissupergrid
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
       z11_to_maggies, cat, maggies, ivarmaggies
       isedfit, isedfit_paramfile, maggies, ivarmaggies, cat.z, $
         result, isedfit_dir=isedfit_dir, isedfit_outfile=isedfit_outfile, $
         supergrid_paramfile=supergrid_paramfile, clobber=clobber, $
         sfhgrid_paramfile=sfhgrid_paramfile, outprefix=outprefix, $
         thissupergrid=thissupergrid, use_redshift=use_redshift
    endif 

; --------------------------------------------------
; make a QAplot
    if keyword_set(qaplot) then begin
       yrange = [31,21.5]
       xrange = [2000,70000]
       xlog = 1
       isedfit_qaplot, isedfit_paramfile, supergrid_paramfile=supergrid_paramfile, $
         isedfit_dir=isedfit_dir, montegrids_dir=montegrids_dir, $
         outprefix=outprefix, clobber=clobber, nsigma=1.1, $
         yrange=yrange, /xlog, thissupergrid=thissupergrid
    endif
    
return
end
