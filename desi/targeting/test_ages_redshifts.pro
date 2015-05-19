
pro test_ages_redshifts, offset=offset, dphase=dphase, ntotal=ntotal, exptime=exptime, outfile=outfile, output_spectra=output_spectra, fluxfac=fluxfac


;; set quicksim params
model_template = 'perfect'
;syfile = '../data/spectra/spec-sky-bright.dat'
syfile = getenv('DESIMODEL')+'data/spectra/spec-sky.dat' ;; use this with Marla+Connie 
desi_fiber_factor = 1 ;; ages = 1.5 arcsec fibers

;; set lunarmodule params
if(not keyword_set(dphase)) then dphase = 14
sepang = 60
ozenang=30
mzenang=30
lunarparam=[dphase,sepang,ozenang,mzenang]

;; number of redshifts
if(not keyword_set(ntotal)) then ntotal = 100

;; exposure time
if(not keyword_set(exptime)) then exptime = 6

;; filenames
if(not keyword_set(outfile)) then outfile = "ages_quicksim_output_nn"

nn = 0
if(keyword_set(offset)) then begin
   nn = offset
   print, 'offset= ', offset
endif

print, 'OUTPUT FILE: ', outfile+strtrim(string(nn),2)+".fits"

;; read in the data
;xh = mrdfits("../ages/primus_ages_parent_v3.0.fits.gz",1)
;x = mrdfits("../ages/primus_ages_parent_v3.0.fits.gz",2)

xh = mrdfits(getenv('DESI_ROOT')+'/spectro/templates/bgs_templates/v1.0/bgs_templates_obs_v1.0.fits', 1)
xs = mrdfits(getenv('DESI_ROOT')+'/spectro/templates/bgs_templates/v1.0/bgs_templates_sdss_v1.0.fits',1)
x = mrdfits(getenv('DESI_ROOT')+'/spectro/templates/bgs_templates/v1.0/bgs_templates_obs_v1.0.fits', 0, hdr)
ngal = n_elements(xh)

;;; remove the aperture correction
;for i=0L, -ngal-1 do begin
;   x.flux[*,i] = x.flux[*,i]/xh[i].apercor
;endfor

;; normalize units to 1E-17 erg/s/cm^2/Ang
x = x/1.0E-17

;; get the output structure
t1 = { id:-1L, mag:0.0, median_sn:0.0, z0:0.0, zfit:0.0, zerr:0.0, rms_flux:0.0 } 
targ = replicate(t1,ngal)

;; get the original redshift
targ.z0 = xh.z

; read the eigentemplates
tt = mrdfits('/Users/ioannis/repos/svn/idlspec2d/templates/spEigenGal-boss_v0_new1.fits',0,thdr)
twave = 10D^(sxpar(thdr,'COEFF0')+dindgen(sxpar(thdr,'NAXIS1'))*sxpar(thdr,'COEFF1'))

for i=0L, ngal-1 do begin
   
   targ[i].id = i

;  flux = x[*,i]*xh[i].infiber_r
   flux = x[*,i]*0.1 ;; TESTING THIS!
   lamz = sxpar(hdr,'CRVAL1')+dindgen(sxpar(hdr,'NAXIS1'))*sxpar(hdr,'CDELT1')

   ;lamz = lamz*(1+targ[i].z0)
   ;flux = flux/(1+targ[i].z0)

   ;; arbitrarily reduce flux by x10
   if(keyword_set(fluxfac)) then begin
      flux = flux/fluxfac
      print, 'FLUXFAC: ',fluxfac
   endif


   ;; get stuff for redshift fitting
   ;wave = lamz
   ;n = n_elements(wave)
   ;dloglam = alog10(wave[n-1]/wave[0])/(n-1)
   ;loglam = findgen(n)*dloglam + alog10(wave[0])
   ;linlam = 10^loglam
   ;flux1= spline(lamz, flux, linlam)
   ;ivar = flux*0 + 0.01
   ;z1 = zfind_tinker(flux1, ivar, xwave=loglam, eigenfile="spEigenGal-53054.fits",zmin=0,zmax=1)
   ;print, z1.z, z1.z_err, targ[i].z0
   ;continue

   if( i eq 0 ) then begin
      desi_quicksim, lamz, flux, exptime=exptime*60, simdat=simdat, model=model_template, skyfile=skyfile, lunarparam=lunarparam, /clear
   endif else begin
      desi_quicksim, lamz, flux, exptime=exptime*60, simdat=simdat, model=model_template, skyfile=skyfile, lunarparam=lunarparam
   endelse
   simdat_original = simdat
   err = randomn(seed,n_elements(simdat.invvar))/(sqrt(simdat.invvar))
;  err = randomn(123,n_elements(simdat.invvar))/(sqrt(simdat.invvar))

   simdat[0].flux = simdat[0].flux + err
   w = where(simdat[0].wavelength ge lamz[0] and simdat[0].wavelength le lamz[n_elements(lamz)-1] and simdat[0].invvar gt 0, n)
   print, total((simdat.flux[w]-err[w])^2)/n, total((simdat.flux[w])^2)/n, total(flux^2)/n
   if(n eq 0) then begin
      targ[i].median_sn = -1
      targ[i].zfit = -1
      targ[i].zerr = -1
      continue
   endif
   targ[i].median_sn = median(simdat[0].snvec[w])

   ;; lets get the S/N in the r-band
   wx = where(simdat.wavelength gt 5600 and simdat.wavelength lt 6700, n1)
   mean_sn = total(simdat.flux[w])/n1

   ;; get stuff for redshift fitting
   wave = simdat[0].wavelength[w]
   n = n_elements(wave)
   dloglam = alog10(wave[n-1]/wave[0])/(n-1)
   ;dloglam = 1.0E-4
   ;n = alog10(wave[n-1]/wave[0])/dloglam+1
   print, 'DLOGLAM', dloglam
   loglam = findgen(n)*dloglam + alog10(wave[0])
   linlam = 10^loglam
   flux = spline(simdat[0].wavelength[w], simdat[0].flux[w], linlam)
      oflux = spline(simdat[0].wavelength[w], simdat_original[0].flux[w], linlam)
   ivar = spline(simdat[0].wavelength[w], simdat[0].invvar[w], linlam)
   print, total((simdat.flux[w]-err[w])^2)/n, total((simdat.flux[w])^2)/n, total(flux^2)/n

   ;; what is size of each pixel in new vector-- scale the variance by
   ;;                                            that value
   dx = dloglam*alog(10)*linlam
   ivar = ivar*(dx/0.5)

   ;; we should print out each spectrum in an ascii file for later
   ;; plotting
   if(keyword_set(output_spectra)) then begin
      openw, lun, "temp_ages_spectra/sp_gal"+padinteger(i,4)+"_exp"+padinteger(exptime,3)+"_dphase10.dat", /get_lun
                                ; make a header for the redshift information
      ;;printf, lun, targ[k].z0, targ[k].zfit, targ[k].zerr
      for i1=0L, n_elements(flux)-1 do begin
         printf, lun, linlam[i1], flux[i1], ivar[i1], oflux[i1]
      endfor
      free_lun, lun
   continue
   endif

   ;;z1 = zfind_tinker(flux, ivar, xwave=loglam, eigenfile="spEigenGal-53054.fits",zmin=0,zmax=3)
;  z1 = zfind(flux, ivar, xwave=loglam, eigenfile="spEigenGal-boss_v0_new1.fits ",zmin=0,zmax=2,zguess1=0.2)
   z1 = zfind_tinker(flux, ivar, xwave=loglam, eigenfile="spEigenGal-boss_v0_new1.fits ",zmin=0,zmax=2,zguess1=0.2)
   targ[i].zfit = z1.z
   targ[i].zerr = z1.z_err

;   ff = flux2maggie(flux, waveimg=10^loglam, mag=mag)
;   flux = x[*,i]*xh[i].infiber_r
;   lamz = sxpar(hdr,'CRVAL1')+dindgen(sxpar(hdr,'NAXIS1'))*sxpar(hdr,'CDELT1')
;   ff = flux2maggie(flux, waveimg=lamz, mag=mag2)
;   print, targ[i].z0, z1.z, z1.z_err, mag[2], mag2[2], xh[i].decam_r, $
;          22.5 - 2.5*alog10(xs[i].modelflux[2]), mean_sn, format='("REDSHIFT ",8(F10.5))'

   djs_plot, 10D^loglam, flux, xsty=3, ysty=3
   djs_oplot, twave*(1+z1.z), tt#z1.theta[0:3], color='orange'



   stop
   cc = get_kbrd(1)

endfor

mwrfits, targ, outfile+strtrim(string(nn),2)+".fits", /create

stop



end
