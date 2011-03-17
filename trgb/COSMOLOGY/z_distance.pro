pro z_distance
;+
; NAME:
;	Z_DISTANCE
;
; PURPOSE:
;	Generate a plot of peculiar velocity versus cz from Peebles'
;	Figure 1 (1988 Apj, 332, 17).
;
; INPUTS:
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;
; PROCEDURES USED:
;	COLORTABLE1
;
; MODIFICATION HISTORY:
;	John Moustakas, 2000 June 27, UCB
;-

        !x.ticklen=0.03
        p_thick = !p.thick & chars = !p.charsize
        charth = !p.charthick & xth = !x.thick & yth = !y.thick
        !p.thick = 2 & !p.charsize=1.8 & !p.charthick=2 & !x.thick=2 & !y.thick=2

	colortable1
        device, decompose=0

	template = {name:	'', $
                    galaxy:	strarr(1), $
                    vc:		fltarr(1), $
                    cz:		fltarr(1)}
        peebles = replicate(template,13)

        peebles[0].galaxy = 'NGC247'
        peebles[0].vc = 157.
        peebles[0].cz = 237.

        peebles[1].galaxy = 'NGC253'
        peebles[1].vc = 244.
        peebles[1].cz = 256.

        peebles[2].galaxy = 'NGC1560'
        peebles[2].vc = 231.
        peebles[2].cz = 215.

        peebles[3].galaxy = 'NGC2403'
        peebles[3].vc = 340.
        peebles[3].cz = 260.

        peebles[4].galaxy = 'NGC3031'
        peebles[4].vc = 169.
        peebles[4].cz = 300.

        peebles[5].galaxy = 'NGC3621'
        peebles[5].vc = 452.
        peebles[5].cz = 505.

        peebles[6].galaxy = 'NGC4144'
        peebles[6].vc = 368.
        peebles[6].cz = 373.

        peebles[7].galaxy = 'NGC4236'
        peebles[7].vc = 212.
        peebles[7].cz = 302.

        peebles[8].galaxy = 'NGC4244'
        peebles[8].vc = 306.
        peebles[8].cz = 342.

        peebles[9].galaxy = 'NGC4826'
        peebles[9].vc = 399.
        peebles[9].cz = 335.

        peebles[10].galaxy = 'NGC5585'
        peebles[10].vc = 488.
        peebles[10].cz = 559.

        peebles[11].galaxy = 'NGC7090'
        peebles[11].vc = 718.
        peebles[11].cz = 541.

        peebles[12].galaxy = 'NGC7793'
        peebles[12].vc = 203
        peebles[12].cz = 306.

; plot the Peebles data

        window, 0, xs=550, ys=550
        plotsym, 0, 1, /fill
        plot, [0,600], [0,600], /nodata, color=1, xsty=3, ysty=3, $
          xtitle='cz km s'+textoidl('^{-1}'), $
          ytitle='v'+textoidl('_{c}')+' km s'+textoidl('^{-1}'), $
          xminor=3, yminor=3, title='Peebles 1988'
        oplot, peebles.cz, peebles.vc, ps=8, color=5
        oplot, [0,600], [0,600], line=0, color=2, thick=2

; read in our results

;       keck_object = sread('/deep1/ioannis/trgb/keck_object.dat')
;       keckdist = keck_object[where(keck_object.distance ne 0.)].distance

;        hst_object = sread('/deep1/ioannis/trgb/hst_object.dat')
;        good = where(hst_object.distance ne 0.)
;        hstdist = hst_object[good].distance
;        hstcz = hst_object[good].v_helio
;
;        h_0 = 69.
;        hstvc = hstcz - h_0*hstdist
;
;        oplot, hstcz, hstvc, ps=2, color=3
        
        stop
        
        !x.ticklen=0
        !p.thick = p_thick & !p.charsize=chars & !p.charthick=charth
        !x.thick = xth & !y.thick=yth

stop
return
end
