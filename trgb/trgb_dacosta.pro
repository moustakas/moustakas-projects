pro trgb_dacosta, gcrgb
;+
; NAME:
;	TRGB_DACOSTA
;
; PURPOSE:
;	Write the I-band magnitudes and V-I colors of the globular cluster
;	RGB's described in Dacosta & Armandroff 1990 to an array of
;	structures. 
;
; INPUTS:
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;	gcrgb	: array of structures for each globular cluster
;		.cluster: globular cluster name
;		.imag	: I-band magnitude
;		.color	; V-I color
;
; COMMON BLOCKS:
;
; SIDE EFFECTS:
;
; RESTRICTIONS:
;
; PROCEDURE:
;
; MODIFICATION HISTORY:
;	John Moustakas, 26 June 2000, UCB
;-

        !x.ticklen=0.03
        p_thick = !p.thick & chars = !p.charsize
        charth = !p.charthick & xth = !x.thick & yth = !y.thick
        !p.thick = 2 & !p.charsize=1.8 & !p.charthick=2 & !x.thick=2 & !y.thick=2

	colortable1

        template = {name:	' ', $
                    cluster:	strarr(1), $
                    imag:	fltarr(28), $
                    color:	fltarr(28)}
        gcrgb = replicate(template,6)

        gcrgb[0].cluster = 'M15'
        gcrgb[0].imag = -[0.224,0.427,0.625,0.817,$
                        1.004,1.186,1.362,1.533,$
                        1.753,1.963,2.163,2.354,$
                        2.535,2.707,2.869,3.021,$
                        3.164,3.329,3.508,3.664,$
                        3.840,3.947,4.020,4.095]
        gcrgb[0].color = [0.866,0.882,0.899,0.915,$
                        0.932,0.948,0.965,0.981,$
                        1.003,1.025,1.047,1.069,$
                        1.091,1.113,1.135,1.157,$
                        1.179,1.209,1.240,1.273,$
                        1.317,1.350,1.377,1.416]

        gcrgb[1].cluster = '47 Tuc'
        gcrgb[1].imag = -[0.131,0.717,1.220,1.651,$
                        2.018,2.329,2.592,2.814,$
                        3.001,3.160,3.294,3.410,$
                        3.510,3.597,3.675,3.745,$
                        3.809,3.866,3.905,3.940,$
                        3.971,3.998,4.019,4.033,$
                        4.039,4.038,4.032,4.021]
        gcrgb[1].color = [0.991,1.060,1.127,1.195,$
                        1.263,1.331,1.399,1.467,$
                        1.535,1.603,1.671,1.739,$
                        1.807,1.875,1.943,2.011,$
                        2.079,2.147,2.198,2.249,$
                        2.300,2.351,2.402,2.453,$
                        2.504,2.538,2.572,2.606]
        
        lincoef = linfit(gcrgb[0].color,gcrgb[0].imag)
        polycoef = poly_fit(gcrgb[0].color,gcrgb[0].imag,2)

        range = minmax(gcrgb[0].color)
        xaxis = findgen((range[1]-range[0])*10^4)/10^4.+range[0]

        linefit = poly(xaxis,lincoef)
        polyfit = poly(xaxis,polycoef)

;        window, 0, xs=550, ys=550
;        plot, [-1.,4.], [-5,0], /nodata, color=1, xsty=3, ysty=3
;        oplot, gcrgb[0].color, gcrgb[0].imag, color=3, line=0
;        oplot, [!x.crange[0],!x.crange[1]], [-4,-4], line=2, color=2, thick=3
;        oplot, xaxis, polyfit, color=5, line=1
;        oplot, xaxis, linefit, color=7, line=2, thick=2
        
; luminosity function
        
        minmag = min(polyfit)
        maxmag = max(polyfit)

        hist = histogram(polyfit,bin=0.1)
        nhist = n_elements(hist)
        xnorm = max(findgen(nhist))/(maxmag-minmag) ; magnitude array
        xmag = findgen(nhist)/xnorm+minmag

;        window, 2, xs=550, ys=550
;        plot, xmag, alog10(hist), xsty=3, ysty=3, color=3, $
;          xtit='I', ytit='Log Number Density'

;stop

        !x.ticklen=0
        !p.thick = p_thick & !p.charsize=chars & !p.charthick=charth
        !x.thick = xth & !y.thick=yth

return
end

