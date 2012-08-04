pro megaspitzer_killerplot, clobber=clobber
; jm12jul29siena - build the killer SED plot; see also
; megaspitzer_killerplot_simulations
    
    path = clash_path(/megaspitzer)
    catfile = path+'killerplot.fits'
    
    cat = read_megaspitzer()
    megaspitzer_to_maggies, cat, refmaggies, refivar
    refmag = maggies2mag(refmaggies,ivar=refivar,magerr=refmagerr)
    ngal = n_elements(cat)
    
    filt = megaspitzer_filterlist()
    nfilt = n_elements(filt)
    isf160 = (where(strmatch(filt,'*f160w*',/fold)))[0]
    isirac = where(strmatch(filt,'*irac*',/fold))
        
; restore the model fits
    isedpath = clash_path(/megaspitzer)+'isedfit/'
    isedfit_sfhgrid_dir = clash_path(/megaspitzer)+'montegrids/'
    sfhgrid_paramfile = getenv('CLASH_DIR')+'/megaspitzer/megaspitzer_sfhgrid.par'

; many the best-fit models are too young to make this a good test, so
; reconstruct the median model, which is older/has a stronger Balmer
; break 
    ndraw = isedfit_ndraw()
    outprefix = 'alldropouts_noirac'
    paramfile = isedpath+'megaspitzer_supergrid01_isedfit.par'
    bestmodel = isedfit_restore(paramfile,ised,iopath=isedpath,$
      isedfit_sfhgrid_dir=isedfit_sfhgrid_dir,/flam,outprefix=outprefix)
    postmodel = bestmodel
    
;    mstar = isedfit_reconstruct_posterior(paramfile,post=post,$
;      isedfit_sfhgrid_dir=isedfit_sfhgrid_dir,iopath=isedpath,$
;      outprefix=outprefix,age=age,Z=Z,tau=tau,sfr0=sfr0,sfrage=sfrage,$
;      chunkindx=chunkindx,modelindx=modelindx,indxage=ageindx,$
;      bigsfr0=bigsfr,bigmass=bigmass,bigsfrage=bigsfrage)
;
;    for ii = 0, ngal-1 do begin
;       get_element, sfrage[*,ii], ised[ii].sfrage_50, this
;       temp = {zobj: ised[ii].zobj, chi2: 0.0, chunkindx: chunkindx[this], $
;         modelindx: modelindx[this], ageindx: ageindx[this], $
;         scale: post[ii].scale[this]}
;       postmodel1 = isedfit_restore(paramfile,in_isedfit=temp,/flam,$
;         isedfit_sfhgrid_dir=isedfit_sfhgrid_dir,iopath=isedpath)
;       if ii eq 0 then postmodel = postmodel1 else postmodel = [postmodel,postmodel1]
;    endfor
    
; build the simulated SEDs    
    sim = replicate({z: 0.0, snr: 0.0, $
      maggies: fltarr(nfilt), ivarmaggies: fltarr(nfilt)},ngal)
    sim.z = cat.z
    sim.snr = 10.0

; assign the S/N of the IRAC bands according to the simulation, then
; add random noise; finally, rescale the IRAC uncertainties so that
; the chosen SNR is preserved
    for ii = 0, ngal-1 do begin
;      bestmaggies = ised[ii].bestmaggies ; best-fit photometry
       bestmaggies = reform(k_project_filters(k_lambda_to_edges(postmodel[ii].wave),$
         postmodel[ii].flux>1D-30,filterlist=filt))

       ivar = ised[ii].ivarmaggies
       ivar[isirac] = (sim[ii].snr/bestmaggies[isirac])^2

       sim[ii].maggies = bestmaggies + randomn(seed,nfilt)/sqrt(ivar)
       sim[ii].ivarmaggies = ivar
       rescale = sim[ii].snr/(sim[ii].maggies[isirac]*sqrt(sim[ii].ivarmaggies[isirac]))
    
       sim[ii].ivarmaggies[isirac] = sim[ii].ivarmaggies[isirac]*rescale^2
;      niceprint, filt, -2.5*alog10(sim[ii].maggies), sim[ii].ivarmaggies, sim[ii].maggies*sqrt(sim[ii].ivarmaggies)
    endfor
; write out the catalog and SEDs
       im_mwrfits, sim, catfile, /clobber, /nogzip

return
end
