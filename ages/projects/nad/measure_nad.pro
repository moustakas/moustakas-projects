
function measure_nad, fiber, restwl, restflux, invvar, continuum, $
      noplot=noplot, title=title, ok = ok, nad_fit = yfit_clip
 
;-----------------------------------------------------------------------------
; Create output structure
;-----------------------------------------------------------------------------

ns = {nad_sn: 0.0, $                ; S/N in region around Na I D
      ew_sub: 0.0, $                ; Index difference between model & data
      ew: 0.0, ew_err: -1.0, $      ; Equivalent Width from integration (A)
      ew_fit: 0.0, $                ; Equivalent Width from model fit (A)
      voff: 0.0, voff_err: -1.0, $  ; Velocity offset (km/s)
      b: 0.0, b_err: -1.0, $        ; Line width (km/s)  
      tau0: 0.0, tau0_err: -1.0, $  ; Optical depth at line center (km/s)
      cf: 0.0, cf_err: -1.0}        ; Covering factor

;------------------------------------------------------------------------------
; Normalize spectrum
;------------------------------------------------------------------------------

norm =  double(restflux / continuum)
norm_err = 1.D / sqrt(invvar) / continuum

; Region fitting is done over
fit_low = 5878.0  ; WATCH OUT FOR He I!!!!
fit_hi = 5905.0 
ok=where(invvar ne 0 and restwl ge fit_low and restwl le fit_hi, nok)
fit_low = min(restwl[ok])
fit_hi = max(restwl[ok])

if nok lt 10 then begin
   print, 'ABORTING FIT!'
   return, ns
end

; Get rid of slight zero-pt offset
cok = where((restwl gt 5820 and restwl lt 5850) or $
            (restwl gt 5920 and restwl lt 5950)) 
corfactor = median(norm[cok])
norm = norm / corfactor
norm_err = norm_err / corfactor
new_cont = continuum * corfactor

;-----------------------------------------------------------------------------
; Get rid of He I 5875.67 line!

emflux = restflux - new_cont

eok = where(restwl gt 5863 and restwl lt 5877)  ; blue side only 
efit_clip = gaussfit(restwl[eok], emflux[eok], ecoef, nterms=3, $
            measure_errors=1/sqrt(invvar[eok]))

gauss_funct, restwl, ecoef, efit 

norm = double(restflux / (new_cont + efit))

;------------------------------------------------------------------------------
; Define starting parameters of Na I D fit
;------------------------------------------------------------------------------

cspeed = 2.99792e5
line_ctr = double([5895.924, 5889.951])

functargs = {dr: [2.0], $ ; dobulet ratio 
    inst_res: [1.D] * fiber.HeI_5876_inst_res, $ ; instrumental resolution
    line_ctr: [line_ctr]} ; line centers of doublet

parinfo = replicate({value:0.D, fixed:0, limited:[0,0], tied:'', $
                     limits:[0.D,0.D]}, 8)

;---------------------------------
; Fit one Na D doublet component

parinfo[0].value = line_ctr[0] ; Starting line center
parinfo[0].limited = [1,1]
parinfo[0].limits = line_ctr[0] * (1 + [-2000.0, 500.0]/cspeed) 

parinfo[1].value = 150.0 ; Starting line b (km/s)
parinfo[1].limited = [1,1]
parinfo[1].limits = [0.65, 800.0]

parinfo[2].value = 1.0 ; Starting tau0
parinfo[2].limited = [1,1]
parinfo[2].limits = [0.1, 5.0]

parinfo[3].value = 0.5 ; Covering factor
parinfo[3].limited = [1,1]
parinfo[3].limits = [0.05, 1.0]

;------------------------------------------------------------------------------
; Fit Interstellar Na D
;------------------------------------------------------------------------------

lfit = mpfitfun('doublet', double(restwl[ok]), norm[ok], norm_err[ok], $
                 functargs=functargs, parinfo=parinfo, perror=perror, $
                 yfit=yfit_clip, maxiter = 300, niter=niter, status=status, $
                 quiet = 1)

print, 'DOUBLETFIT ITERATIONS: ', strtrim(niter, 2)
print, 'DOUBLETFIT EXIT STATUS: ', strtrim(status, 2)

;------------------------------------------------------------------------------
; Populate output structure
;------------------------------------------------------------------------------

; S/N 
ns.nad_sn = median(norm[ok]/norm_err[ok])

; Calculate EW in 3 different ways
ns.ew_sub = fiber.na_d_abs - fiber.na_d_abs_model 

ns.ew = integral(restwl[ok], 1 - norm[ok], fit_low, fit_hi)
ns.ew_err = sqrt(integral(restwl[ok], norm_err[ok]^2, fit_low, fit_hi)) 

ns.ew_fit = integral(restwl[ok], 1 - yfit_clip, fit_low, fit_hi)

; Velocity offset (km/s)
ns.voff = (lfit[0] - line_ctr[0]) / line_ctr[0] * cspeed
ns.voff_err = perror[0] / line_ctr[0] * cspeed

; Intrinsic line width (km/s)
ns.b = lfit[1]
ns.b_err = perror[1]

; Optical depth a line center
ns.tau0 = lfit[2]
ns.tau0_err = perror[2]

; Covering factor
ns.cf = lfit[3]
ns.cf_err = perror[3]

print, ns.ew, ns.voff, format='(3F8.2)'

;------------------------------------------------------------------------------
; Diagnostic Plot
;------------------------------------------------------------------------------

if not keyword_set(noplot) then begin

  pagemaker, nx=1, ny=2, pos=pos_top, xpage=10, ypage=7.5, /norm, $
             xmargin = 0.5, ymargin=0.5
  pagemaker, nx=2, ny=2, pos=pos_bot, xpage=10, ypage=7.5, /norm, $
             xmargin = 0.5, ymargin=0.5

  ; Continuum fit to full spectrum
  plot, restwl, restflux, /xs, yr=[0, max(continuum) * 1.2], /ys, $
        xr = [3600, 7200], title=title, pos=pos_top[*,0]

  djs_oplot, restwl, continuum, color='red', thick=2
  xyouts, 0.10, 0.90, 'S/N = ' + string(fiber.sn_median, format='(I3)'), /norm

  ; Region around Na D 
  plot, restwl, restflux, xr = [5820, 5950], yr=[0.7, 1.25], /xs, /ys, $
        thick=2, pos=pos_bot[*,2], /noerase
  
  djs_oplot, restwl, new_cont, color='red', thick=2

  ycfit = yfit_clip*new_cont[ok]
  smooth_model = smooth(medsmooth(new_cont, 100), 30, /edge_truncate)
  yfit = yfit_clip*smooth_model[ok]
  djs_oplot, restwl[ok], ycfit, color='green', thick=2 
  djs_oplot, restwl[ok], yfit, color='blue', thick=2 

  ;eok = where(restwl gt 5850 and restwl lt 5885)
  djs_oplot, restwl[eok], efit[eok] + new_cont[eok], color='cyan', thick=2

  ; Normalized interstellar line (show both components)
  plot, restwl, norm, /xs, /ys, yr = [0.85, 1.05], $
        xr = [5840, 5950], thick=2, pos=pos_bot[*,3], /noerase

  oplot, [5000, 6000], [1, 1]  
  djs_oplot, restwl[ok], yfit_clip, color='green', thick=2, linestyle=0

  xyouts, 0.55, 0.16, 'V_off = ' + string(ns.voff, form='(F7.1)'), /norm
  xyouts, 0.55, 0.12, 'Na D EW = ' + string(ns.ew_fit, form='(F6.1)'), /norm

endif

return, ns

end
