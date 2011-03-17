FUNCTION fit_kurucz_models, flux, ivar, filterlist=filterlist
   
      ;;fstarfile = '$PRIMUS_DIR/etc/primus_atlas9_stds.fits'
    fstarfile = '$IDLSPEC2D_DIR/etc/kurucz_stds_v5.fit'

    fstarstr = mrdfits(fstarfile, 1)
    fstarspec = mrdfits(fstarfile, 0, fstarhdr)
      
      ;;fstarwave = mrdfits(fstarfile, 2)
    npix = n_elements(fstarspec[*,0])
    fstarwave = 10d^(sxpar(fstarhdr, 'crval1') + $
                       sxpar(fstarhdr, 'cd1_1')*findgen(npix))
   
   nstar = n_elements(fstarstr)
   nphot = n_elements(flux[*,0])
   ngal = n_elements(flux[0,*])
   fstar_flux = dblarr(nphot, nstar)
   wave = k_lambda_eff(filterlist=strtrim(filterlist,2))
   
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;;
   ;;  READ the fstar photometry
   ;; 
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   FOR istar = 0L, n_elements(fstarstr)-1L DO BEGIN
      fstar_flux[*,istar] = k_project_filters(k_lambda_to_edges(fstarwave),$
        reform(fstarspec[*,istar]),filterlist=strtrim(filterlist,2))
   ENDFOR

   chi2 = fltarr(nstar, ngal)+1E6
   scaled_flux = fltarr(nphot, nstar, ngal)
   scale = fltarr(nstar, ngal)
   
   ;;;;;;;;;;;;;;;;;;;;;;
   ;;
   ;; Now, find the best fit star
   ;;
   ;;;;;;;;;;;;;;;;;;;;;
   
   FOR istar = 0, nstar - 1 DO BEGIN
      FOR igal = 0, ngal -1 DO BEGIN
; require photometry in at least three bandpasses
         good = where((flux[*,igal] gt 0.0) and (ivar[*,igal] gt 0.0),ngood)
         if (ngood ge 3) then begin
            scale[istar,igal] = $
              total(fstar_flux[*,istar]*flux[*,igal]*ivar[*,igal], $
              /double) / $
              total(fstar_flux[*,istar]*fstar_flux[*,istar]*ivar[*,igal], /double)
            scaled_flux[*,istar,igal] = fstar_flux[*,istar]*scale[istar, igal]
            chi2[istar, igal] = total((flux[*,igal]-scaled_flux[*,istar,igal])^2 $
              *ivar[*,igal])
         endif
      ENDFOR
   ENDFOR

   lambda = fstarwave
   spec = fltarr(n_elements(lambda), ngal)
   chi2min = fltarr(ngal)
   best = intarr(ngal)
   FOR igal = 0, ngal -1 DO BEGIN
      chi2min[igal] = min(chi2[*,igal], kmin)
      spec[*,igal] = fstarspec[*,kmin]*scale[kmin, igal]
      best[igal] = kmin
   ENDFOR

   nfilt = n_elements(filterlist)
   output = {$
     flag:          0, $ ; 0=good, 1=bad
     kurucz_model: '', $
     kurucz_feh:   -999.0, $
     kurucz_teff:  -999.0, $
     kurucz_g:     -999.0, $
     primus_maggies: fltarr(nfilt), $
     primus_maggiesivar: fltarr(nfilt), $
     kurucz_filters: strarr(nfilt), $
     kurucz_maggies:  fltarr(nfilt), $
     kurucz_chi2min: 0.0, $
     best: -1, $
     lambda: float(lambda), $
     spec: lambda*0.0}
   output = replicate(output, ngal)

   output.primus_maggies = flux
   output.primus_maggiesivar = ivar
   output.kurucz_filters = filterlist
   output.kurucz_chi2min = chi2min

   flag = where(output.kurucz_chi2min eq 1E6,nflag,$
     comp=good,ncomp=ngood)
   if (nflag ne 0) then output[flag].flag = 1

   if (ngood ne 0) then begin
      output[good].spec = spec[*,good]
      output[good].best = best[good]
      
      output[good].kurucz_model = fstarstr[best[good]].model
      output[good].kurucz_feh = fstarstr[best[good]].feh
      output[good].kurucz_teff = fstarstr[best[good]].teff
      output[good].kurucz_g = fstarstr[best[good]].g
      
      for jj = 0, ngood-1 do output[good[jj]].kurucz_maggies = $
        scaled_flux[*,best[good[jj]],good[jj]]
   endif

return, output
END

