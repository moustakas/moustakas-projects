;+
; NAME:
;   TRGB_MAXLIKE
; PURPOSE:
;   to run likelihood analysis on a starlist, looking for TRGB in I
;   band.  Data based either on keck run or on HST study.
; CALLING SEQUENCE:
;   TRGB_MAXLIKE, objname, [minmag = minmag, maxmag = maxmag, $
;                  alpha = alpha, /hst, /halo, /core, /ccut]
; INPUTS:
;   objname -- name of galaxy
; KEYWORD PARAMETERS:
;   alpha -- set this to the slope of the RGB (default: 0.3)
;   hst -- select if hst object, otherwise keck
;   halo -- select if desire to read halo of object only
;   core -- select if desire to read core of object only
;   ccut -- select to read the color-cut photometry of selected regions
;   minmag -- min magnitude (I=20 default)
;   maxmag -- max magnitude (default is max of histogram binned at 0.05 mag)
; PROCEDURES USED:
;   ERR_FUNC
;   TRGB_LFMODEL
; MODIFICATION HISTORY:
;  jm, md  23jun00
;  bjm 29sep00
;-
pro TRGB_MAXLIKE, objname, minmag = minmag, maxmag = maxmag, $
                  alpha = alpha, $
                  hst = hst, halo = halo, core = core, ccut = ccut, ps = ps

; read in starlist ----------------------------------------------------------- 
TRGB_READATA, objname, datapath, data, info, $
              hst = hst, halo = halo, core = core, ccut = ccut

  imags = data[3, *]
  PRINT, N_ELEMENTS(imags), ' stars input for '+info.truename & PRINT

; Set Magnitude range --------------------------------------------------------
  IF NOT KEYWORD_SET(minmag) THEN  minmag = 20.0 

  hist = HISTOGRAM(imags, bin = 0.05, min = minmag)
  m = FINDGEN(N_ELEMENTS(hist))*0.05 + minmag
  hmax = MAX(hist, hm)
  mmax = m(hm)

  IF NOT KEYWORD_SET(maxmag) THEN  maxmag = mmax

; Calculate the likelihood --------------------------------------------------- 
bmag = ROUND((imags - minmag)*100.) ;binned magnitudes in range minmag, maxmag
numberbins = 100*(maxmag-minmag)
bmag = bmag[WHERE((bmag ge 0) and (bmag lt numberbins))] ; subset of stars in selected magnitude range
 
Ntotal = N_ELEMENTS(bmag)
lnlike = FLTARR(50,20,6)
PRINT, Ntotal,' stars within range', minmag, ' to ', maxmag

; Produce the error function to boxcar smooth each model with
       mmag = FINDGEN(500)*0.01 + minmag
       error = ERR_FUNC(objname, mmag, hst = hst, halo = halo, core = core)
       nn = 500
       dm = 100
       ewidth = ROUND(error*dm*1.44) ;cube root of 3 to get correct width sigma
       eplus = (INDGEN(nn) + ewidth) < (nn-1)
       a1 = WHERE(eplus EQ nn-1)
       sa = SIZE(a1)
       s1 = sa[1]
       eminus = (INDGEN(nn) - ewidth) > 0
       eminus[a1[0]:a1[s1-1]] = eminus[a1[0]]
       nsmooth = eplus - eminus + 1

FOR l = 0, 5 DO BEGIN           ;loop over beta
   beta = l/10. +.5             ;slope of the bright end of LF
   FOR k = 0, 19 DO BEGIN       ;loop over discontinuity strength
      cutfactor = k/20. +.1
      FOR i = 0, 49 DO BEGIN    ;loop over cutoff magnitude
         mc = i/10.             ; magnitude cutoff in range minmag, minmag+5
         gg = TRGB_LFMODEL_bm(error, mc, cutfactor, eminus, eplus, nsmooth, a1,$
                           s1, alpha = alpha, beta = beta)
         aloggg = ALOG(gg)
         lnlike[i, k, l] = -Ntotal*ALOG(TOTAL(gg[0:numberbins-1])) + $
                            TOTAL(aloggg(bmag))
      ENDFOR     
   ENDFOR
ENDFOR

lnlike = lnlike - MAX(lnlike) ; normalize the likelihood

; regenerate best fit ---------------------------------------------------------

mag = FINDGEN(50)/10. + minmag
step = FINDGEN(20)/20. + 0.1
bstep = FINDGEN(6)/10. + 0.5

dummy = MAX(lnlike, i) ; N.B. lnlike has dimension (50,20,6)
ibeta = i/(50*20)
betaf = ibeta/10. +.5
cut = FIX((i-ibeta*(50*20))/50.)
cutf = cut/20. + 0.1
magcut = (i mod 50)/10. + minmag  
PRINT, 'best fit paramters: m_cut = ',magcut, ' cutf= ',cutf, ' beta',betaf

gg = TRGB_LFMODEL(error, magcut-minmag, cutf, eminus, eplus, nsmooth, a1, s1,$
                  alpha = alpha, beta = betaf)  

GETELEMENT_VECTOR, m, magcut, x1
x2 = ROUND(100*(magcut-minmag))
best_fit = gg*hist(x1)/gg(x2)

; Plots ----------------------------------------------------------------------
COLORTABLE1
WINDOW, 0, xs = 500, ys = 700

PLOT, m, (hist > 1), psym = 10, xrange = [minmag, mmax], /ylog, $
  title = info.truename, ytitle = 'log number', $
  position = [0.1,0.5,0.95,0.95]
OPLOT, mmag, best_fit, lines = 3, color = 7, thick = 2
OPLOT, [minmag, minmag], [1, !y.crange[1]], $
     line = 2, thick = 2, color = 2
OPLOT, [maxmag, maxmag], [1, !y.crange[1]], $
     line = 2, thick = 2, color = 2
   IF (KEYWORD_SET(halo) and KEYWORD_SET(ccut)) THEN $
legend, ['TRGB = '+STRN(magcut,format='(F6.2)'), $
         'cutfactor = '+STRN(cutf, format='(F4.2)'), $
         'brightend slope = '+STRN(betaf, format = '(F4.2)'), 'halo,ccut'], $
        charsize = 1.2, textcolor=[7,7,7, 7], box=0, /clear $
   ELSE $
   IF (KEYWORD_SET(core) and KEYWORD_SET(ccut)) THEN $
legend, ['TRGB = '+STRN(magcut,format='(F6.2)'), $
         'cutfactor = '+STRN(cutf, format='(F4.2)'), $
         'brightend slope = '+STRN(betaf, format = '(F4.2)'), 'core,ccut'], $
        charsize = 1.2, textcolor=[7,7,7, 7], box=0, /clear $
   ELSE $
   IF KEYWORD_SET(ccut) THEN $
legend, ['TRGB = '+STRN(magcut,format='(F6.2)'), $
         'cutfactor = '+STRN(cutf, format='(F4.2)'), $
         'brightend slope = '+STRN(betaf, format = '(F4.2)'), 'color-cut'], $
        charsize = 1.2, textcolor=[7,7,7, 7], box=0, /clear $
   ELSE $
   IF KEYWORD_SET(halo) THEN $
legend, ['TRGB = '+STRN(magcut,format='(F6.2)'), $
         'cutfactor = '+STRN(cutf, format='(F4.2)'), $
         'brightend slope = '+STRN(betaf, format = '(F4.2)'), 'halo'], $
        charsize = 1.2, textcolor=[7,7,7, 7], box=0, /clear $
   ELSE $
   IF KEYWORD_SET(core) THEN $
legend, ['TRGB = '+STRN(magcut,format='(F6.2)'), $
         'cutfactor = '+STRN(cutf, format='(F4.2)'), $
         'brightend slope = '+STRN(betaf, format = '(F4.2)'), 'core'], $
        charsize = 1.2, textcolor=[7,7,7, 7], box=0, /clear $
   ELSE $
legend, ['TRGB = '+STRN(magcut,format='(F6.2)'), $
         'cutfactor = '+STRN(cutf, format='(F4.2)'), $
         'brightend slope = '+STRN(betaf, format = '(F4.2)')], $
        charsize = 1.2, textcolor=[7,7,7], box=0, /clear 

lnlikecf = lnlike[*,*,ibeta]
CONTOUR, lnlikecf, mag, step, levels=[-3,-2,-1], c_colors = [1, 7, 2], $
  ytitle = 'log cut amplitude', xrange = [minmag, mmax], $
  position = [0.1,0.3,0.95,0.5], /NOERASE

lnlikeb = REFORM(lnlike[*,cut,*])
CONTOUR, lnlikeb, mag, bstep, levels = [-3,-2,-1], c_colors = [1, 7, 2], $
  xtitle = 'magnitude', ytitle='slope of bright end', $
  position = [0.1,0.1,0.95,0.3], xrange = [minmag, mmax], /NOERASE


IF KEYWORD_SET(ps) THEN BEGIN
PS_OPEN, 'ml', /portrait

PLOT, m, (hist > 1), psym = 10, xrange = [minmag, mmax], /ylog, $
  title = info.truename, ytitle = 'log number', $
  position = [0.1,0.5,0.95,0.95]
OPLOT, mmag, best_fit, lines = 3, color = 7, thick = 2
OPLOT, [minmag, minmag], [1, !y.crange[1]], $
     line = 2, thick = 2, color = 2
OPLOT, [maxmag, maxmag], [1, !y.crange[1]], $
     line = 2, thick = 2, color = 2
   IF (KEYWORD_SET(halo) and KEYWORD_SET(ccut)) THEN $
legend, ['TRGB = '+STRN(magcut,format='(F6.2)'), $
         'cutfactor = '+STRN(cutf, format='(F4.2)'), $
         'brightend slope = '+STRN(betaf, format = '(F4.2)'), 'halo,ccut'], $
        charsize = 1.2, textcolor=[7,7,7, 7], box=0, /clear $
   ELSE $
   IF (KEYWORD_SET(core) and KEYWORD_SET(ccut)) THEN $
legend, ['TRGB = '+STRN(magcut,format='(F6.2)'), $
         'cutfactor = '+STRN(cutf, format='(F4.2)'), $
         'brightend slope = '+STRN(betaf, format = '(F4.2)'), 'core,ccut'], $
        charsize = 1.2, textcolor=[7,7,7, 7], box=0, /clear $
   ELSE $
   IF KEYWORD_SET(ccut) THEN $
legend, ['TRGB = '+STRN(magcut,format='(F6.2)'), $
         'cutfactor = '+STRN(cutf, format='(F4.2)'), $
         'brightend slope = '+STRN(betaf, format = '(F4.2)'), 'color-cut'], $
        charsize = 1.2, textcolor=[7,7,7, 7], box=0, /clear $
   ELSE $
   IF KEYWORD_SET(halo) THEN $
legend, ['TRGB = '+STRN(magcut,format='(F6.2)'), $
         'cutfactor = '+STRN(cutf, format='(F4.2)'), $
         'brightend slope = '+STRN(betaf, format = '(F4.2)'), 'halo'], $
        charsize = 1.2, textcolor=[7,7,7, 7], box=0, /clear $
   ELSE $
   IF KEYWORD_SET(core) THEN $
legend, ['TRGB = '+STRN(magcut,format='(F6.2)'), $
         'cutfactor = '+STRN(cutf, format='(F4.2)'), $
         'brightend slope = '+STRN(betaf, format = '(F4.2)'), 'core'], $
        charsize = 1.2, textcolor=[7,7,7, 7], box=0, /clear $
   ELSE $
legend, ['TRGB = '+STRN(magcut,format='(F6.2)'), $
         'cutfactor = '+STRN(cutf, format='(F4.2)'), $
         'brightend slope = '+STRN(betaf, format = '(F4.2)')], $
        charsize = 1.2, textcolor=[7,7,7], box=0, /clear 

lnlikecf = lnlike[*,*,ibeta]
CONTOUR, lnlikecf, mag, step, levels=[-3,-2,-1], c_colors = [1, 7, 2], $
  ytitle = 'log cut amplitude', xrange = [minmag, mmax], $
  position = [0.1,0.3,0.95,0.5], /NOERASE

lnlikeb = REFORM(lnlike[*,cut,*])
CONTOUR, lnlikeb, mag, bstep, levels = [-3,-2,-1], c_colors = [1, 7, 2], $
  xtitle = 'magnitude', ytitle='slope of bright end', $
  position = [0.1,0.1,0.95,0.3], xrange = [minmag, mmax], /NOERASE

PS_CLOSE
ENDIF


  return
end
