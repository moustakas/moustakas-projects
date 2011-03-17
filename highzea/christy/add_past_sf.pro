

function add_past_sf, time, sfr, model, eval_age, wlrange = wlrange, $
         wave = wave, mass = mass, ok = ok, ssp_nlyc=ssp_nlyc, $
         total_nlyc = total_nlyc

; model.age [nage] in yrs!
; model.flux [npix, nage] in L_sun/A for 1 M_sun of stars!
; model.wave [npix]

; to avoid loss of resolution at old ages, rebin 10 x


; Compute SFR at models of each age
pop_age = eval_age - time
linterp, pop_age, sfr, model.age, isfr

;plot, pop_age, sfr, psym=4, syms=1, xr=[0, eval_age]
;oplot, model.age, isfr, color=!red

; compensate for resolution loss
nsfr = isfr
nage = n_elements(model.age)

for ii = 1, nage - 2 do begin
  dtl = (model.age[ii] - model.age[ii - 1])/2
  dth = (model.age[ii + 1] - model.age[ii])/ 2
  bin = where(pop_age gt model.age[ii] - dtl and $
              pop_age le model.age[ii] + dth)
  if bin[0] eq -1 then continue
  nsfr[ii] = avg(sfr[bin])
endfor
;oplot, model.age, nsfr, color=!blue
;wait, 2

if not keyword_set(wlrange) then $
  wlrange = [min(model.wave), max(model.wave)]

ok = where(model.wave ge wlrange[0] and $
           model.wave le wlrange[1], npix)
wave = model.wave[ok]
spec = fltarr(npix)

dt0 = min(model.age) ; first time step

; Compute integral over each pixel
for ii = 0, npix - 1 do begin
  ;spec[ii] = tsum(model.age, isfr*model.flux[ok[ii],*], 0, bad[0]-1)
   spec[ii] = im_integral(model.age, isfr*model.flux[ok[ii],*], $
              min(model.age), eval_age)
   ; add in fist fraction of a Myr not covered by integral 
   spec[ii] = spec[ii] + isfr[0]*model.flux[ok[ii],0]*dt0
endfor

mass = im_integral(model.age, isfr, min(model.age), eval_age)
mass = mass + isfr[0]*dt0

if keyword_set(ssp_nlyc) then begin
  total_nlyc = im_integral(double(model.age), isfr*ssp_nlyc, $
                           min(model.age), eval_age)
  total_nlyc = total_nlyc + isfr[0] * ssp_nlyc[0]*dt0
endif

return, spec

end
