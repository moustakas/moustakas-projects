
pro mass_oh, mass, oh, figname = figname, medm = medm, medoh = medoh, $
    sigoh = sigoh, ohfit = ohfit, ohother = ohother, nohist = nohist, $
    show = show

;******************************************************************************
; The Mass-Metallicity relation
;******************************************************************************

if not keyword_set(figname) then figname = 'mass_oh'

;-----------------------------------
; Plot 2d histogram

plot2dh, mass, oh, [7.7, 11.7], [8.0,9.5], figname, img = img, pixthresh=5, $
;  ytitle = '12 + log(O/H)', xtitle = 'log(M!L*!N/M!3!DO!N!X)', /ys, /xs, $
  ytitle = '12 + log(O/H)', xtitle = 'log(M!L*!N)', /ys, /xs, $
  charsize = 1.5, charthick =2, thick = 3, xthick=3, ythick=3
;xyouts, 9.996, 7.858, '.', charthick=3, charsize=1  ; dot in M_sun

medxbin, mass, oh, 0.1, medx=medm, medy=medoh, sigy=sigoh, thresh=100, distr=d
;djs_oploterr, medm, medoh, yerr=sigoh, psym=8, thick = 4
oplot, medm, medoh, psym=8, thick=5
;oplot, medm, medoh, psym=8, symsize = 1.2, thick=5
for i = 0, 3 do oplot, medm, d[i,*], thick=5

if keyword_set(show) then $
  oplot, mass[show], oh[show], psym=4, syms=0.3, color=!red

;-----------------------
; Create table

mz_table, medm, medoh, d, xvar = '$\log$(M$_{*}/$M$_{\sun}$)', $
          title = 'The Mass--Metallicity Relation', tabname='tb3.tex'


;-----------------------------------
; Fit M-Z relation with Legendre polynomial

;coef = robust_poly_fit(medm, medoh, 2, ohfit)
coef = mpfitfun('polyleg', medm, medoh, sigoh, [1, 1, 1], yfit=ohfit, /quiet)
oplot, medm, ohfit, color=!red, thick=5

mfid = findgen(60) /10 + 6
ohmod = polyleg(mfid, coef)

cnew = [coef[0] - coef[2]/2, coef[1], coef[2] * 3./2]
print, 'M-Z Coeff: ', cnew

;ohmod = -1.501 + 1.843 * mfid - 0.05340 * (3.0 * mfid^2 - 1) / 2.0
ohmod = -1.492 + 1.847 * mfid - 0.08026 * mfid^2 
;oplot, mfid, ohmod, color=!green

a = value_to_index(medm, 10.0)
print, 'Scatter in M-Z relation at M=10^10:', sigoh[a]
print, 'Median scatter in M-Z relation:', median(d[2,*] - d[1,*])/2, $
       median(d[3,*] - d[0,*])/2


;-------------------------------------
; Jarle's broken powerlaw fit

;coef2 = mpfitfun('brokenpowerlaw', medm, medoh, sigoh, $
;                [2.0, 5.0, -2.0, 10.0, 2.0], yfit=yfitbp, quiet=1)
;print, coef2
;oplot, mfid, brokenpowerlaw(mfid, coef2), color=!green
;oplot, mfid, polyleg(mfid, coef), color=!red

;------------------------------------------------------
; If desired overplot fit to different measure of O/H

if keyword_set(ohother) then begin
  
  ok = where(ohother ne 0)

  medxbin, mass[ok], ohother[ok], 0.1, medx=medmo, medy = medoho, $
           sigy = sigoho, thresh=100
  ocoef = mpfitfun('polyleg', medmo, medoho, sigoho, [1, 1, 1], $
         yfit=ohfito, /quiet)
  ;print, ocoef
  ;oplot, medmo, medoho, psym=8, thick=5, color=!dpink
  ;oplot, medmo, medoho, psym=8, symsize = 1.2, thick=5
  oplot, medmo, ohfito, lines = 2, thick=5
endif

;--------------------------------------
; Plot residuals

ohfit = polyleg(mass, coef)
if not keyword_set(nohist) then begin
  nplothist, oh - ohfit, bin=0.01, xhist, yhist, /noplot
  plot, xhist, yhist, /noerase, position=[0.7, 0.18, 0.90, 0.42], psym=10, $
      xtitle='12 + log(O/H) Residuals', charsize = 0.8, $
      ytitle='Number of Galaxies', thick=3, /norm, xr=[-1, 1], /xs

  xyouts, 0.1, 2200, '!9s!X = ' + string(robust_sigma(oh - ohfit, /zero), $
        format = '(F5.2)')

  ; --------------------
  ; Brute force calculation of residual sigma
  ;resid = oh - ohfit
  ;resid = resid[where(finite(resid))]
  ;resid = resid[sort(resid)]
  ;nx = n_elements(resid)
  ;linterp, findgen(nx)/nx, resid, [0.025, 0.16, 0.5, 0.84, 0.975], residp
  ;print, residp, format='(5F6.2)'

endif


stopprint

end

