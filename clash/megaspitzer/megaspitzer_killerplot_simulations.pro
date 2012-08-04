pro megaspitzer_killerplot_simulations, clobber=clobber
; jm12jul29siena - build the killer SED plot

; the simulations produced by this piece of code are great, but we
; ended up not using them in the proposal
    
    path = clash_path(/megaspitzer)
    cat = read_megaspitzer()
    megaspitzer_to_maggies, cat, refmaggies, refivar
    refmag = maggies2mag(refmaggies,ivar=refivar,magerr=refmagerr)
    
    filt = megaspitzer_filterlist()
    nfilt = n_elements(filt)
    isf160 = (where(strmatch(filt,'*f160w*',/fold)))[0]
    isirac = where(strmatch(filt,'*irac*',/fold))
        
;   ssp = im_read_fsps(/flam)
    ssp = im_read_fsps(met='Z0.0039',/flam)
    nwave = n_elements(ssp.wave)

    igmgrid = mrdfits(getenv('IMPRO_DIR')+'/etc/igmtau_grid.fits.gz',1)

    h100 = 0.7D
    pc10 = 3.085678D19 ; fiducial distance [10 pc in cm]
    light = 2.99792458D18       ; speed of light [A/s]

    catfile = path+'killerplot_simulations.fits'
    if im_file_test(catfile,clobber=clobber) then return
    
    nmod = 500
    zmin = 4.0 & zmax = 9.5
    taumin = 0D & taumax = 0.5D
    snrmin = 0.1 & snrmax = 30.0      ; for [ch1], [ch2]
    minmgal = 1D8 & maxmgal = 1D11    ; total galaxy mass
    minage = 0.1D & maxage = 1.2D
;   minage = getage(15D) & maxage = 1.2D
    
    sim = replicate({z: 0.0, tau: 0D, age: 0.0, mgal: 0.0, $
      mass: 0.0, snr: 0.0, $
      maggies: fltarr(nfilt), ivarmaggies: fltarr(nfilt), $
      wave: fltarr(nwave), flux: fltarr(nwave)},nmod)
    sim.z = randomu(seed,nmod)*(zmax-zmin)+zmin
    sim.snr = 10.0
;   sim.snr = 10.0^(randomu(seed,nmod)*(alog10(snrmax)-alog10(snrmin))+alog10(snrmin))
    sim.tau = randomu(seed,nmod)*(taumax-taumin)+taumin
    sim.age = randomu(seed,nmod)*(maxage-minage)+minage
    sim.mgal = alog10(randomu(seed,nmod)*(maxmgal-minmgal)+minmgal)

    dlum = pc10*10D^(lf_distmod(sim.z,omega0=0.3D,$ ; [cm]
      omegal0=0.7D)/5D)/h100
    
; convolve with the SFH, redshift, add IGM, then errors
    for ii = 0, nmod-1 do begin

       sfh = {tau: sim[ii].tau, delayed: 1, nburst: 0, tautrunc: 0, maxage: maxage}
;      sfh1 = isedfit_reconstruct_sfh(sfh,/deb)

       csp1 = isedfit_convolve_sfh(ssp,mstar=ssp.mstar,info=sfh,$
         cspmstar=cspmstar,time=sim[ii].age)
       flux1 = 10.0^sim[ii].mgal*csp1*(pc10/dlum[ii])^2D/(1+sim[ii].z)

       sim[ii].mass = sim[ii].mgal+alog10(cspmstar)
       sim[ii].wave = ssp.wave*(1+sim[ii].z)

       windx = findex(igmgrid.wave,sim[ii].wave)
       zindx = findex(igmgrid.zgrid,sim[ii].z)
       igm = interpolate(igmgrid.igm,windx,zindx,/grid,missing=1.0)
       sim[ii].flux = flux1*igm
       
; convolve with the HST filters; to add uncertainties, choose randomly
; from the Megaspitzer dropouts catalog
       maggies1 = reform(k_project_filters(k_lambda_to_edges(sim[ii].wave),$
         sim[ii].flux>1D-30,filterlist=filt))
       mag1 = -2.5*alog10(maggies1)
;      niceprint, mag1, filt

;      djs_plot, sim[ii].wave, -2.5*alog10((sim[ii].flux*sim[ii].wave^2/light)>1D-50)-48.6, $
;        /xlog, xr=[1500,1D5];, yr=[35,20]
;      djs_plot, sim[ii].wave, flux1, /xlog, xr=[1500,1D5], /ylog
;      djs_oplot, sim[ii].wave, sim[ii].flux, color='orange'

       get_element, refmag[isf160,*], mag1[isf160], this
       ivar1 = refivar[*,this]
       ivar = ivar1

; assign the S/N of the IRAC bands according to the simulation, then
; add random noise; finally, rescale the IRAC uncertainties so that
; the chosen SNR is preserved
       ivar[isirac] = (sim[ii].snr/maggies1[isirac])^2
       maggies = maggies1 + randomn(seed,nfilt)/sqrt(ivar)
       
       minerr = replicate(0.1,nfilt)
       minerr[isirac] = 0.0 ; note!
       k_minerror, maggies, ivar, minerr

       sim[ii].maggies = maggies
       sim[ii].ivarmaggies = ivar
       rescale = sim[ii].snr/(sim[ii].maggies[isirac]*sqrt(sim[ii].ivarmaggies[isirac]))
       
       sim[ii].ivarmaggies[isirac] = sim[ii].ivarmaggies[isirac]*rescale^2
;      niceprint, filt, -2.5*alog10(sim[ii].maggies), sim[ii].ivarmaggies, sim[ii].maggies*sqrt(sim[ii].ivarmaggies)
    endfor

; write out the catalog and SEDs
    im_mwrfits, sim, catfile, /clobber, /nogzip
    
return
end
