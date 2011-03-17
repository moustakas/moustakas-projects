;+
; NAME:
;   sdss2ages
; PURPOSE:
;   take SDSS data and return NDWFS/AGES photometry at some redshift
; CALLING SEQUENCE:
;   agesphot= sdss2ages(sdss_redshift, ages_redshift, [, nmgy=, ivar=, mag=, err=, $
;                  calibobj=, tsobj=, flux=, chi2=, rmaggies=, $
;                  omaggies=, vname=, oivar=, mass=, mtol= ]
; INPUTS:
;   sdss_redshift - [N] redshifts of input 
;   ages_redshift - [N] redshifts of desired output
;   calibobj - [N] photoop-style structure, containing:
;                  .PETROFLUX[5]
;                  .PETROFLUX_IVAR[5]
;                  .MODELFLUX[5]
;                  .MODELFLUX_IVAR[5]
;                  .PSFFLUX[5]
;                  .PSFFLUX_IVAR[5]
;                  .EXTINCTION[5]
;   tsobj - [N] opdb-style structure, containing:
;                  .PETROCOUNTS[5]
;                  .PETROCOUNTSERR[5]
;                  .COUNTS_MODEL[5]
;                  .COUNTS_MODELERR[5]
;                  .PSFCOUNTS[5]
;                  .PSFCOUNTSERR[5]
;                  .REDDENING[5]
;   nmgy, ivar - [5, N] nanomaggies, Galactic-reddening corrected, and inverse
;                variance of same
;   mag, err - [5, N] asinh magnitudes, Galactic-reddening corrected and
;              errors of same
; OPTIONAL INPUTS:
;   flux - use this version of the fluxes ('PETRO', 'MODEL', or 'PSF')
;          [defaults to 'PETRO'] if tsobj or calibobj keywords are
;          used 
;   vname - name of fit to use (defaults to 'default')
; OUTPUTS:
;   bwrik - [4, N] apparent magnitudes in BRI (AB)
;   mtol - [5, N] mass-to-light ratios from model in each band
;   mass - [N] total mass from model in each band
; OPTIONAL OUTPUTS:
;   coeffs - coefficients of fit
;   chi2 - chi^2 of fit
;   rmaggies - [5, N] reconstructed maggies from the fit (ugriz)
;   omaggies, oivar - [5, N] maggies and inverse variances used for fit
;                           (after extinction, AB correction, etc)  (ugriz)
; COMMENTS:
;   This is a simple wrapper on kcorrect.pro. It keeps a version of
;   rmatrix and zvals in memory to save time.
;
;   You must specify nmgy,ivar OR mag,err OR calibobj OR tsobj. If
;   nmgy or mag, make sure they are AB calibrated and Galactic
;   extinction corrected.
;
;   Uses sdss_to_maggies to convert tsobj or calibobj structure to
;   AB, Galactic extinction corrected maggies. Passes optional
;   argument "flux" to sdss_to_maggies.
;
;   For v4_0b templates and later, coefficients are in units of: 
; 
;     1 solar mass / (D/10pc)^2 
;
;   That is, sum the coefficients and multiply by (D/10pc)^2 to get
;   TOTAL INTEGRATED STAR FORMATION. (In fact, for Omega0=0.3 and
;   OmegaL0=0.7, this is what the "mass" keyword returns). Note that
;   the total integrated star formation DIFFERS from the current
;   stellar mass --- which is returned in the mass and mtol variables.
; REVISION HISTORY:
;   07-Apr-2005  Mike Blanton, NYU
;-
;------------------------------------------------------------------------------

function sdss2ages, sdss_redshift, ages_redshift, nmgy=nmgy, ivar=ivar, $
  mag=mag, err=err, calibobj=calibobj, tsobj=tsobj, flux=flux, chi2=chi2, $
  coeffs=coeffs, rmaggies=rmaggies, omaggies=omaggies, oivar=oivar, vname=vname, $
  mass=mass, mtol=mtol, agesphot_ivar=agesphot_ivar, silent=silent, q0=in_q0, $
  q1=in_q1, qz0=in_qz0, out_filterlist=out_filterlist, in_filterlist=in_filterlist

    common com_sdss2ages, rmatrix, zvals, band_shift

    if (n_params() lt 1 OR $
      (((keyword_set(nmgy) eq 0 OR keyword_set(ivar) eq 0)) AND $
      ((keyword_set(mag) eq 0 OR keyword_set(err) eq 0)) AND $
      (n_tags(calibobj) eq 0) AND $
      (n_tags(tsobj) eq 0))) $
    then begin
       doc_library, 'sdss2ages'
       return, -1
    endif 

    if (n_elements(in_q0) eq 0) then in_q0=0.
    if (n_elements(in_q1) eq 0) then in_q1=0.
    if (n_elements(in_qz0) eq 0) then in_qz0=0.

    q0 = in_q0
    q1 = in_q1
    qz0 = in_qz0

    if (n_elements(in_filterlist) eq 0) then in_filterlist = sdss_filterlist()
    if (n_elements(out_filterlist) eq 0) then out_filterlist = ages_filterlist()

    ngal = n_elements(sdss_redshift)
    nfilt_out = n_elements(out_filterlist)
    nfilt_in = n_elements(in_filterlist)
    
;   t0 = systime(1)
    splog, 'Note: using SDSS_KCORRECT (requires just using the SDSS filters!)'
    kcdum = sdss_kcorrect(sdss_redshift,nmgy=nmgy,ivar=ivar,mag=mag,$
      err=err,calibobj=calibobj,tsobj=tsobj,flux=flux,chi2=chi2,$
      coeffs=coeffs,rmaggies=rmaggies,omaggies=omaggies,oivar=oivar,$
      vname=vname,mass=mass,mtol=mtol,band_shift=band_shift,/silent)
;   print, systime(1)-t0

; calculate the preliminaries
    if (not keyword_set(rmatrix) OR not keyword_set(zvals)) then begin
       if (not keyword_set(vmatrix) OR not keyword_set(lambda)) then $
         k_load_vmatrix, vmatrix, lambda, vfile=vfile, lfile=lfile, $
         vpath=vpath, vname=vname
       k_projection_table, rmatrix, vmatrix, lambda, zvals, out_filterlist, $ 
         zmin=zmin, zmax=zmax, nz=nz, filterpath=filterpath, /silent
    endif

; reconstruct the magnitudes as observed by AGES
    k_reconstruct_maggies,coeffs, ages_redshift, $
      reconstruct_maggies, rmatrix=rmatrix, zvals=zvals, $
      /silent
    obands = lindgen(nfilt_out)#replicate(1L,ngal)

    closest = 1 ; NOTE!
    if keyword_set(closest) then begin
       lambda_in = k_lambda_eff(filterlist=in_filterlist)
       lambda_out = k_lambda_eff(filterlist=out_filterlist)
       for i=0L, ngal-1L do begin
          for j=0L, n_elements(lambda_out)-1L do begin
             dmin=min(abs(lambda_in/(1.+sdss_redshift[i])- $
               lambda_out[j]/(1.0+ages_redshift[i])), imin)
             obands[j, i]= imin
          endfor
       endfor
    endif

    offset = fltarr(nfilt_out,ngal)
    for i = 0L, ngal-1L do $
      for j = 0L, nfilt_out-1L do offset[j,i] = $
      reconstruct_maggies[j,i]/rmaggies[obands[j,i],i]
    offset = 2.5*alog10(offset)

    agesphot = fltarr(nfilt_out,ngal)
    agesphot_ivar = fltarr(nfilt_out,ngal)
    dm_sdss = lf_distmod(sdss_redshift,omega0=omega0,omegal0=omegal0)
    dm_ages = lf_distmod(ages_redshift,omega0=omega0,omegal0=omegal0) - $ ; allow evolution!!
      k_evolve(ages_redshift*0.0,ages_redshift,q0,q1,qz0)

    for j=0L, nfilt_out-1L do $
      agesphot[j,*]=-2.5*alog10(reconstruct_maggies[j,*])-dm_sdss+dm_ages 

    for i = 0L, nfilt_in-1L do begin
       for j = 0L, nfilt_out-1L do begin
          ifrom = where(obands[j,*] eq i, nfrom)
          if (nfrom gt 0) then begin
             ig = where(oivar[i,ifrom] gt 0.0 and $
               omaggies[i,ifrom] gt 0.0,ng)
             if (ng gt 0) then begin
                ig = ifrom[ig]
                agesphot[j,ig] = -2.5*alog10(omaggies[i,ig]) - $
                  dm_sdss[ig]+dm_ages[ig]-offset[j,ig]
                agesphot_ivar[j,ig] = omaggies[i,ig]^2*oivar[i,ig]*(0.4*alog(10))^2
             endif
          endif
       endfor
    endfor

return, agesphot
end
