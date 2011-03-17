pro TRGB_DERIV, objname, binsize=binsize, minmag=minmag, maxmag=maxmag, $
          hst=hst, halo=halo, core=core, ccut=ccut, $
          ps=ps, noplot = noplot, trgb = trgb, slope = slope

        IF NOT KEYWORD_SET(binsize) THEN binsize = 0.03
;        bin size should be in the range 0.005 to 0.05. around 0.03
;        seems to give the best results in tests with bintest.pro

; read in the data ----------------------------------------------------------
        trgb_readata, objname, datapath, data, infobase, $
                      hst = hst, halo = halo, core = core, ccut = ccut

        nstars = n_elements(data[0,*])

        IF NOT KEYWORD_SET(minmag) THEN minmag = min(data[3,*])
        IF NOT KEYWORD_SET(maxmag) THEN maxmag = max(data[3,*])

; generate a binned luminosity function --------------------------------------

        trgb_lfunction, data[3,*], data[4,*], lf, binsize=binsize, $
                        minmag = minmag, maxmag = maxmag

; generate the continuous luminosity function -------------------------------

phi = FLTARR(lf.nhist)
err = FLTARR(lf.nhist)
FOR j = 0L, lf.nhist-1L DO BEGIN
     phi[j] = TOTAL((1./(SQRT(2.*!PI)*lf.merr))*$
                    EXP(-0.5*((lf.mags-lf.mag[j])/lf.merr)^2))
ENDFOR

; Test that the normalization is correct -------------------------------------

Print, 'Integral of Phi(m) =', (TOTAL(phi)*binsize)/TOTAL(lf.hist)
stop
; Define smooth magnitude error function -------------------------------------

sig = ERR_FUNC(objname, lf.mag, hst=hst, halo=halo, core=core)
sig = 3*sig ; a cheat to beat down the noise in the bright end of the lf

; Define the Response Functions ---------------------------------------------
resplin=FLTARR(lf.nhist)
resplog=FLTARR(lf.nhist)
FOR k=1,lf.nhist-2 DO BEGIN
   GETELEMENT_VECTOR, lf.mag, (lf.mag[k]+sig[k]), p1
   GETELEMENT_VECTOR, lf.mag, (lf.mag[k]-sig[k]), p2
   resplin[k] = phi[p1] - phi[p2]                         ; linear
   resplog[k] = ALOG10(phi[p1] > 1) - ALOG10(phi[p2] > 1) ; log
ENDFOR

; Poisson Noise -------------------------------------------------------------

    smooth = SMOOTH2(phi,0.2/binsize); smooth the continous luminosity function

        noise = SQRT(smooth > 1)
        respnoise = noise*SQRT(2)

        weightedlin = resplin/respnoise
        weightedlog = resplog*respnoise

; Find the TRGB -------------------------------------------------------------

dummy = MAX(phi, phimax)
tweightedlin = WHERE((lf.mag LT lf.mag[phimax]))
tweightedlog = WHERE((lf.mag LT lf.mag[phimax]))

IF KEYWORD_SET(trgb) THEN BEGIN
   GETELEMENT_VECTOR, lf.mag, (trgb-0.5), r1
   GETELEMENT_VECTOR, lf.mag, (trgb+0.5), r2
   a1 = WHERE(weightedlin EQ MAX(weightedlin[r1:r2]))
   trgblin = lf.mag[a1] 
   a2 = WHERE(weightedlog EQ MAX(weightedlog[r1:r2]))
   trgblog = lf.mag[a2]
ENDIF ELSE BEGIN
   a1 = WHERE(weightedlin EQ MAX(weightedlin(tweightedlin)))
   trgblin = lf.mag[a1] 
   a2 = WHERE(weightedlog EQ MAX(weightedlog(tweightedlog)))
   trgblog = lf.mag[a2]
ENDELSE

PRINT, 'Linear TRGB = ', trgblin
PRINT, 'Log TRGB = ', trgblog

; Determine RGB power law ---------------------------------------------------- 

        IF NOT KEYWORD_SET(slope) THEN slope = 0.3
        x3 = phimax - FIX(0.5/binsize)
        rgb = phi[x3]*10^(slope*(lf.mag - lf.mag[x3]))

        diff = rgb - phi
        logdiff = ALOG10(rgb > 1) - ALOG10(phi > 1)

        GETELEMENT_VECTOR, lf.mag, (trgblog + 0.5), x4
        rgb2 = phi[x4]*10^(slope*(lf.mag - lf.mag[x4]))

;      if rgb and rgb2 do not overlap then the real trgb has not
;      been detected or there is some strange feature at the
;      completeness limit of the luminosity function


; Plots ----------------------------------------------------------------------
IF NOT KEYWORD_SET(noplot) THEN BEGIN

   COLORTABLE1
   WINDOW, 0, xs = 500, ys = 700

   PLOT, lf.mag, phi, YRANGE = [1, MAX(phi)], $
     xsty = 3, ysty = 3, color = 1, $
     position = [0.1, 0.5, 0.95, 0.95],  $
     TITLE = 'Linear Luminosity Function and Response for '+infobase.truename,$
     YTITLE = 'Number Density'
   OPLOT, lf.mag, rgb, line = 3, thick = 2, color = 5
   OPLOT, lf.mag, rgb2, line = 3, thick = 2, color = 6
   OPLOT, [trgblin, trgblin], [!y.crange[0], !y.crange[1]], $
     line = 2, thick = 2, color = 2
   legend, ['Stars = '+strn(nstars), $
            'Binsize = '+strn(lf.binsize, format = '(F4.2)'), $
            'TRGB = '+strn(trgblin, format = '(F6.2)'), $
            'RGB'], line = [0, 0, 0, 3], colors = [0, 0, 0, 5], $
            charsize = 1.2, textcolor = [1, 1, 1, 1], box = 0, /clear 
   OPLOT, lf.mag, smooth, color = 7

   PLOT, lf.mag, weightedlin, /NoErase, $
    YRANGE = [MIN(weightedlin(tweightedlin)), MAX(weightedlin(tweightedlin))],$
     xsty = 3, ysty = 3, $
     position = [0.1, 0.25, 0.95, 0.5], $
     YTITLE = 'Weighted Response'
   OPLOT, [trgblin, trgblin], [!y.crange[0], !y.crange[1]], $
     line = 2, thick = 2, color = 2

   PLOT, lf.mag, diff, $
     YRANGE = [MIN(diff), MAX(diff[0:phimax])], /NoErase, $
     xsty = 3, ysty = 3, $
     position = [0.1, 0.1, 0.95, 0.25], $
     YTITLE = 'RGB difference', XTITLE = 'I Magnitude'
   OPLOT, [trgblin, trgblin], [!y.crange[0], !y.crange[1]], $
     line = 2, thick = 2, color = 2

   WINDOW, 10, xs = 500, ys = 700

   PLOT, lf.mag, phi, YRANGE = [1, MAX(phi)], /YLOG, $
     xsty = 3, ysty = 3, color = 1, $
     position = [0.1, 0.5, 0.95, 0.95],  $
TITLE = 'Logarithmic Luminosity Function and Response for '+infobase.truename,$
     YTITLE = 'Number Density' 
   OPLOT, lf.mag, rgb, line = 3, thick = 2, color = 5
   OPLOT, lf.mag, rgb2, line = 4, thick = 3, color = 6
   OPLOT, [trgblog, trgblog], [1, phi[a2]], line = 2, thick = 3, color = 2
   OPLOT, lf.mag, smooth, color = 7
   IF (KEYWORD_SET(halo) and KEYWORD_SET(ccut)) THEN $
      legend, ['Stars = '+strn(nstars), $
               'Binsize = '+strn(lf.binsize, format = '(F4.2)'), $
               'TRGB = '+strn(trgblog, format = '(F6.2)'), 'halo,ccut'], $
        charsize = 1.2, textcolor = [1, 1, 1, 1], box = 0, /clear $ 
   ELSE $
   IF (KEYWORD_SET(core) and KEYWORD_SET(ccut)) THEN $
      legend, ['Stars = '+strn(nstars), $
               'Binsize = '+strn(lf.binsize, format = '(F4.2)'), $
               'TRGB = '+strn(trgblog, format = '(F6.2)'), 'core,ccut'], $
        charsize = 1.2, textcolor = [1, 1, 1, 1], box = 0, /clear $
   ELSE $
   IF KEYWORD_SET(ccut) THEN $
     legend, ['Stars = '+strn(nstars), $
              'Binsize = '+strn(lf.binsize, format = '(F4.2)'), $
              'TRGB = '+strn(trgblog, format = '(F6.2)'), 'color-cut'], $
         charsize = 1.2, textcolor = [1, 1, 1, 1], box = 0, /clear $
   ELSE $
   IF KEYWORD_SET(halo) THEN $
      legend, ['Stars = '+strn(nstars), $
               'Binsize = '+strn(lf.binsize, format = '(F4.2)'), $
               'TRGB = '+strn(trgblog, format = '(F6.2)'), 'halo'], $
        charsize = 1.2, textcolor = [1, 1, 1, 1], box = 0, /clear $
   ELSE $
   IF KEYWORD_SET(core) THEN $
      legend, ['Stars = '+strn(nstars), $
               'Binsize = '+strn(lf.binsize, format = '(F4.2)'), $
               'TRGB = '+strn(trgblog, format = '(F6.2)'), 'core'], $
        charsize = 1.2, textcolor = [1, 1, 1, 1], box = 0, /clear $
   ELSE $
   legend, ['Stars = '+strn(nstars), $
            'Binsize = '+strn(lf.binsize, format = '(F4.2)'), $
            'TRGB = '+strn(trgblog, format = '(F6.2)')], $
            charsize = 1.2, textcolor = [1, 1, 1], box = 0, /clear 


   PLOT, lf.mag, weightedlog, /NoErase, $
YRANGE = [MIN(weightedlog(tweightedlog)), MAX(weightedlog(tweightedlog))], $
     xsty = 3, ysty = 3, $
     position = [0.1, 0.25, 0.95, 0.5], $
     YTITLE = 'Weighted Response'
   OPLOT, [trgblog, trgblog], [!y.crange[0], !y.crange[1]], $
     line = 2, thick = 2, color = 2

   PLOT, lf.mag, logdiff, $
     YRANGE = [MIN(logdiff), MAX(logdiff[0:phimax])], /NoErase, $
     xsty = 3, ysty = 3, $
     position = [0.1, 0.1, 0.95, 0.25], $
     YTITLE = 'RGB difference', XTITLE = 'I Magnitude'
   OPLOT, [trgblog, trgblog], [!y.crange[0], !y.crange[1]], $
     line = 2, thick = 2, color = 2

ENDIF

; Postscript output -----------------------------------------------------------
IF KEYWORD_SET(ps) THEN BEGIN

PS_OPEN, 'lf', /portrait
PLOT, lf.mag, phi, YRANGE = [1, MAX(phi)],  PSYM = 10, /YLOG, $
xsty = 3, ysty = 3, color = 1, $
position = [0.1,0.6,0.92,0.92],  $
TITLE = 'Logarithmic Luminosity Function and Response for '+infobase.truename,$
YTITLE = 'Number Density' 
OPLOT, lf.mag, phi, color = 7
OPLOT, lf.mag, rgb, line=3, thick=2, color=5
OPLOT, [trgblog,trgblog], [1,phi[a2]], line=2, thick=3, color=2
   IF (KEYWORD_SET(halo) and KEYWORD_SET(ccut)) THEN $
      legend, ['Stars = '+strn(nstars), $
               'Binsize = '+strn(lf.binsize, format = '(F4.2)'), $
               'TRGB = '+strn(trgblog, format = '(F6.2)'), 'halo,ccut'], $
        charsize = 1.2, textcolor = [1, 1, 1, 1], box = 0, /clear $ 
   ELSE $
   IF (KEYWORD_SET(core) and KEYWORD_SET(ccut)) THEN $
      legend, ['Stars = '+strn(nstars), $
               'Binsize = '+strn(lf.binsize, format = '(F4.2)'), $
               'TRGB = '+strn(trgblog, format = '(F6.2)'), 'core,ccut'], $
        charsize = 1.2, textcolor = [1, 1, 1, 1], box = 0, /clear $
   ELSE $
   IF KEYWORD_SET(ccut) THEN $
     legend, ['Stars = '+strn(nstars), $
              'Binsize = '+strn(lf.binsize, format = '(F4.2)'), $
              'TRGB = '+strn(trgblog, format = '(F6.2)'), 'color-cut'], $
         charsize = 1.2, textcolor = [1, 1, 1, 1], box = 0, /clear $
   ELSE $
   IF KEYWORD_SET(halo) THEN $
      legend, ['Stars = '+strn(nstars), $
               'Binsize = '+strn(lf.binsize, format = '(F4.2)'), $
               'TRGB = '+strn(trgblog, format = '(F6.2)'), 'halo'], $
        charsize = 1.2, textcolor = [1, 1, 1, 1], box = 0, /clear $
   ELSE $
   IF KEYWORD_SET(core) THEN $
      legend, ['Stars = '+strn(nstars), $
               'Binsize = '+strn(lf.binsize, format = '(F4.2)'), $
               'TRGB = '+strn(trgblog, format = '(F6.2)'), 'core'], $
        charsize = 1.2, textcolor = [1, 1, 1, 1], box = 0, /clear $
   ELSE $
   legend, ['Stars = '+strn(nstars), $
            'Binsize = '+strn(lf.binsize, format = '(F4.2)'), $
            'TRGB = '+strn(trgblog, format = '(F6.2)')], $
            charsize = 1.2, textcolor = [1, 1, 1], box = 0, /clear 

PLOT, lf.mag, resplog, $
YRANGE = [MIN(resplog(tweightedlog)), MAX(resplog(tweightedlog))], /NoErase, $
   xsty = 3, ysty = 3, $
   position = [0.1,0.4,0.92,0.6], $
YTITLE = 'Response Function'
OPLOT, [trgblog,trgblog], [!y.crange[0],!y.crange[1]], line=2, thick=2, color=2

PLOT, lf.mag, weightedlog, /NoErase, $
YRANGE = [MIN(weightedlog(tweightedlog)), MAX(weightedlog(tweightedlog))], $
   xsty = 3, ysty = 3, $
   position = [0.1,0.2,0.92,0.4], $
YTITLE = 'Weighted Response'
OPLOT, [trgblog,trgblog], [!y.crange[0],!y.crange[1]], line=2, thick=2, color=2

PLOT, lf.mag, logdiff, $
YRANGE = [MIN(logdiff) , MAX(logdiff[0:phimax])], /NoErase, $
   xsty = 3, ysty = 3, $
   position = [0.1,0.1,0.92,0.2], $
YTITLE = 'RGB difference', XTITLE = 'I Magnitude'
OPLOT, [trgblog,trgblog], [!y.crange[0],!y.crange[1]], line=2, thick=2, color=2

PS_CLOSE

ENDIF


  return
end
