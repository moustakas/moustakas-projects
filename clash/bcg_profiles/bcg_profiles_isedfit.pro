pro bcg_profiles_isedfit, prelim=prelim, models=models, isedfit=isedfit, $
  qaplot=qaplot, clobber=clobber, bcgmodel=bcgmodel
; jm13jun26siena
    
    prefix = 'bcg_profiles'
    isedfit_dir = getenv('CLASH_DATA')+'/bcg_profiles/isedfit/'
    isedfit_paramfile = isedfit_dir+prefix+'_paramfile.par'
    sfhgrid_paramfile = isedfit_dir+prefix+'_sfhgrid.par'
    supergrid_paramfile = isedfit_dir+prefix+'_supergrid.par'

; get the sample (not all of these have BCG profiles yet)
    clash = rsex(getenv('CLASH_DIR')+'/clash_sample.sex')
    ncl = n_elements(clash)

    use_redshift = clash.z ; use a custom redshift array
    zmin = min(use_redshift)
    zmax = max(use_redshift)
    nzz = n_elements(use_redshift)

    filterlist = bcgimf_filterlist(instr=instr)
    wfc3ir = where(instr eq 'wfc3ir')
    nfilt = n_elements(filterlist)
    
; supergrid priors
    sfhgrid = 1
    supergrid = 1
    synthmodels = 'fsps'
    igm = 0
    
    imf = 'salp'
    redcurve = -1
;   redcurve = 2 ; Milky Way
    
; SFHgrid priors    
    nage = 10 ; 20
    nmonte = 200 ; 500
    Z = [0.004,0.04]
    minage = 1.0
    maxage = 11.0
    AV = [0.0,0.0]
    delay_tau = [0.01,3.0]
    
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

;; simple tau
;       write_sfhgrid_paramfile, sfhgrid_paramfile, params, clobber=0, /append, $
;         sfhgrid=2, nage=nage, nmonte=nmonte, Z=Z, minage=minage, $
;         maxage=maxage, AV=AV, tau=[0.01,3.0], flatAV=flatAV, delayed=0
;; simple tau, maximally old
;       write_sfhgrid_paramfile, sfhgrid_paramfile, params, clobber=0, /append, $
;         sfhgrid=3, nage=nage, nmonte=nmonte, Z=Z, minage=10.5, $
;         maxage=10.5, AV=AV, tau=[0.01,3.0], flatAV=flatAV, delayed=0

       write_supergrid_paramfile, supergrid_paramfile, supergrid=supergrid, $
         sfhgrid=sfhgrid, synthmodels=synthmodels, imf=imf, redcurve=redcurve, $
         clobber=clobber

       build_montegrids, sfhgrid_paramfile, supergrid_paramfile=supergrid_paramfile, $
         isedfit_dir=isedfit_dir, clobber=clobber;, thissupergrid=1
    endif

; --------------------------------------------------
; build the models
    if keyword_set(models) then begin
       isedfit_models, isedfit_paramfile, isedfit_dir=isedfit_dir, $
         supergrid_paramfile=supergrid_paramfile, clobber=clobber, $
         use_redshift=use_redshift;, thissupergrid=1
    endif

; --------------------------------------------------
; do the fitting, looping on each cluster
    if keyword_set(isedfit) then begin
       for ic = 0, ncl-1 do begin
          phot = read_bcg_profiles(clash[ic].shortname,these_filters=filterlist)
          if size(phot,/type) eq 8 then begin
             outprefix = strtrim(clash[ic].shortname,2)

; add a 2% floor on the uncertainties and then fit each radial bin
; independently (effectively treating each radial bin as a different
; galaxy); require good photometry in at least three of the WFC3/IR
; filters
             good = where(total(phot[wfc3ir].ivarmaggies gt 0,2) ge 3.0,ngood)
             radius = phot[wfc3ir[0]].sma[good] ; [kpc]
             maggies = transpose(phot.maggies[good])
             ivarmaggies = transpose(phot.ivarmaggies[good])
             k_minerror, maggies, ivarmaggies, replicate(0.02,nfilt)
             
             isedfit, isedfit_paramfile, maggies, ivarmaggies, replicate(clash[ic].z,ngood), $
               result, isedfit_dir=isedfit_dir, isedfit_outfile=isedfit_outfile, $
               supergrid_paramfile=supergrid_paramfile, clobber=clobber, $
               sfhgrid_paramfile=sfhgrid_paramfile, outprefix=outprefix

; rewrite the results with the radius appended
             result = mrdfits(isedfit_outfile+'.gz',1,/silent)
             result = struct_addtags(result,replicate({radius: 0.0},ngood))
             result.radius = radius
             im_mwrfits, result, isedfit_outfile, /clobber

             splog, outprefix, alog10(total(10^result.mass_50))
stop             
          endif 
       endfor
    endif

; --------------------------------------------------
; make a QAplot
    if keyword_set(qaplot) then begin
       isedfit_qaplot, isedfit_paramfile, supergrid_paramfile=supergrid_paramfile, $
         isedfit_dir=isedfit_dir, montegrids_dir=montegrids_dir, $
         galaxy=galaxy1, outprefix=outprefix, clobber=clobber, /xlog;, thissupergrid=1
    endif
    
    
return
end
    
