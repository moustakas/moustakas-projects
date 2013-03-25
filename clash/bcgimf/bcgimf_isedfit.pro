pro bcgimf_isedfit, prelim=prelim, models=models, isedfit=isedfit, $
  qaplot=qaplot, clobber=clobber, bcgmodel=bcgmodel
; jm13mar15siena
    
    prefix = 'bcgimf'
    isedfit_dir = getenv('BCGIMF_DATA')+'/'
    isedfit_paramfile = isedfit_dir+prefix+'_paramfile.par'
    sfhgrid_paramfile = isedfit_dir+prefix+'_sfhgrid.par'
    supergrid_paramfile = isedfit_dir+prefix+'_supergrid.par'

; global parameters    
    filterlist = bcgimf_filterlist()
    
    zbcg = 0.2233 ; Coe+12
    zmin = 0.22
    zmax = 0.23
    nzz = 3

; SFHgrid priors    
    nage = 20
    nmonte = 500
    Z = [0.004,0.04]
    minage = 1.0
    maxage = 11.0
    AV = [0.0,0.0]

; supergrid parameters    
    sfhgrid = [1,2,2,3]
    supergrid = [1,2,3,4]
    synthmodels = ['fsps','fsps','bc03','fsps']

    imf = 'salp'
    redcurve = -1
;   redcurve = 2 ; Milky Way
    
; --------------------------------------------------
; do the preliminaries: build the parameter files and the Monte Carlo
; grids
    if keyword_set(prelim) then begin
       write_isedfit_paramfile, filterlist, prefix=prefix, minz=zmin, $
         maxz=zmax, nzz=nzz, igm=0, isedfit_dir=isedfit_dir, clobber=clobber

; delayed
       write_sfhgrid_paramfile, sfhgrid_paramfile, params, clobber=clobber, $
         sfhgrid=1, nage=nage, nmonte=nmonte, Z=Z, minage=minage, $
         maxage=maxage, AV=AV, tau=[0.01,3.0], flatAV=flatAV, delayed=1
; simple tau
       write_sfhgrid_paramfile, sfhgrid_paramfile, params, clobber=0, /append, $
         sfhgrid=2, nage=nage, nmonte=nmonte, Z=Z, minage=minage, $
         maxage=maxage, AV=AV, tau=[0.01,3.0], flatAV=flatAV, delayed=0
; simple tau, maximally old
       write_sfhgrid_paramfile, sfhgrid_paramfile, params, clobber=0, /append, $
         sfhgrid=3, nage=nage, nmonte=nmonte, Z=Z, minage=10.5, $
         maxage=10.5, AV=AV, tau=[0.01,3.0], flatAV=flatAV, delayed=0

       write_supergrid_paramfile, supergrid_paramfile, supergrid=supergrid, $
         sfhgrid=sfhgrid, synthmodels=synthmodels, imf=imf, redcurve=redcurve, $
         clobber=clobber

       build_montegrids, sfhgrid_paramfile, supergrid_paramfile=supergrid_paramfile, $
         isedfit_dir=isedfit_dir, clobber=clobber, thissupergrid=4
    endif

; --------------------------------------------------
; build the models
    if keyword_set(models) then begin
       isedfit_models, isedfit_paramfile, isedfit_dir=isedfit_dir, $
         supergrid_paramfile=supergrid_paramfile, clobber=clobber, thissupergrid=4
    endif

; --------------------------------------------------
; do the fitting
    if keyword_set(isedfit) then begin
       if keyword_set(bcgmodel) then begin
          phot = mrdfits(isedfit_dir+'bcgimf_phot_bcgmodel.fits.gz',1) 
          outprefix = 'bcgmodel'
       endif else begin
          phot = mrdfits(isedfit_dir+'bcgimf_phot.fits.gz',1)
       endelse
       nradius = n_elements(phot)

       isedfit, isedfit_paramfile, phot.maggies, phot.ivarmaggies, $
         replicate(zbcg,nradius), result, isedfit_dir=isedfit_dir, $
         supergrid_paramfile=supergrid_paramfile, clobber=clobber, $
         sfhgrid_paramfile=sfhgrid_paramfile, outprefix=outprefix, thissupergrid=4
    endif

; --------------------------------------------------
; make a QAplot
    if keyword_set(qaplot) then begin
       isedfit_qaplot, isedfit_paramfile, supergrid_paramfile=supergrid_paramfile, $
         isedfit_dir=isedfit_dir, montegrids_dir=montegrids_dir, $
         galaxy=galaxy1, outprefix=outprefix, clobber=clobber, /xlog, thissupergrid=4
    endif
    
    
return
end
    
