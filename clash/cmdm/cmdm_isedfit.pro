pro cmdm_isedfit, prelim=prelim, models=models, isedfit=isedfit, $
  qaplot=qaplot, clobber=clobber, fordoron=fordoron
; jm13jul18siena
    
    prefix = 'cmdm'
    datapath = getenv('IM_PROJECTS_DIR')+'/clash/cmdm/'
    isedfit_dir = datapath+'isedfit/'
    isedfit_paramfile = isedfit_dir+prefix+'_paramfile.par'
    sfhgrid_paramfile = isedfit_dir+prefix+'_sfhgrid.par'
    supergrid_paramfile = isedfit_dir+prefix+'_supergrid.par'

; read the catalog and then specify the redshifts and filters
    cat = mrdfits(datapath+'cmdm_cat.fits.gz',1)
    ngal = n_elements(cat)

    filterlist = clash_subaru_filterlist()
    nfilt = n_elements(filterlist)

    zmin = 0.40 ; focused on MACS1206 for now
    zmax = 0.48
    nzz = 10

; supergrid priors
    sfhgrid = 1
    supergrid = 1
    synthmodels = 'fsps'
    igm = 0
    
    imf = 'salp'
    redcurve = 1 ; Charlot & Fall
    
; SFHgrid priors    
    nage = 20
    nmonte = 500
    Z = [0.004,0.04]
    minage = 0.1
    maxage = 9.2
    delay_tau = [0.001,5.0]

; --------------------------------------------------
; do the preliminaries: build the parameter files and the Monte Carlo
; grids
    if keyword_set(prelim) then begin
       write_isedfit_paramfile, filterlist, prefix=prefix, minz=zmin, $
         maxz=zmax, nzz=nzz, igm=igm, isedfit_dir=isedfit_dir, clobber=clobber

; delayed
       write_sfhgrid_paramfile, sfhgrid_paramfile, params, clobber=clobber, $
         sfhgrid=1, nage=nage, nmonte=nmonte, Z=Z, minage=minage, $
         maxage=maxage, tau=delay_tau, /delayed

       write_supergrid_paramfile, supergrid_paramfile, supergrid=supergrid, $
         sfhgrid=sfhgrid, synthmodels=synthmodels, imf=imf, redcurve=redcurve, $
         clobber=clobber

       build_montegrids, sfhgrid_paramfile, supergrid_paramfile=supergrid_paramfile, $
         isedfit_dir=isedfit_dir, clobber=clobber, thissupergrid=1
    endif

; --------------------------------------------------
; build the models
    if keyword_set(models) then begin
       isedfit_models, isedfit_paramfile, isedfit_dir=isedfit_dir, $
         supergrid_paramfile=supergrid_paramfile, clobber=clobber
    endif

; --------------------------------------------------
; do the fitting!
    if keyword_set(isedfit) then begin
       clash_subaru_to_maggies, cat, maggies, ivarmaggies
          
       isedfit, isedfit_paramfile, maggies, ivarmaggies, cat.zobj, $
         result, isedfit_dir=isedfit_dir, isedfit_outfile=isedfit_outfile, $
         supergrid_paramfile=supergrid_paramfile, clobber=clobber, $
         sfhgrid_paramfile=sfhgrid_paramfile, outprefix=outprefix
    endif 

; --------------------------------------------------
; make a QAplot
    if keyword_set(qaplot) then begin
       indx = shuffle_indx(ngal,num=50)
       isedfit_qaplot, isedfit_paramfile, supergrid_paramfile=supergrid_paramfile, $
         isedfit_dir=isedfit_dir, montegrids_dir=montegrids_dir, index=indx, $
         galaxy=galaxy1, outprefix=outprefix, clobber=clobber, /xlog
    endif

; --------------------------------------------------
; package the results for Doron
    if keyword_set(fordoron) then begin
       ised = read_isedfit(isedfit_paramfile,isedfit_dir=isedfit_dir,$
         supergrid_paramfile=supergrid_paramfile)
       ngal = n_elements(ised)
       mstar = isedfit_reconstruct_posterior(isedfit_paramfile,$
         supergrid_paramfile=supergrid_paramfile,thissupergrid=thissupergrid,$
         isedfit_dir=isedfit_dir,montegrids_dir=montegrids_dir)

       ndraw = isedfit_ndraw()
       out = struct_addtags(cat,struct_trimtags(ised,select=['mass_50','mass_err']))
       out = struct_addtags(out,replicate({pofm: fltarr(ndraw)},ngal))
       out.pofm = mstar

       wsex, struct_trimtags(out,except='pofm'), outfile=datapath+'cmdm_isedfit.txt'
       openw, lun, datapath+'cmdm_isedfit_pofm.txt', /get_lun
       for ii = 0L, ngal-1 do printf, lun, 
       free_lun, lun

stop       
       
       im_mwrfits, out, datapath+'cmdm_isedfit.fits', clobber=clobber

       
       
stop       
    endif
       
    
return
end
    
