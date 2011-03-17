pro plot_smf, objname, lf, phi, resplin, resplog, trgblin, trgblog, infobase, nstars, ps=ps

	!x.ticklen=0.02
        p_thick = !p.thick & chars = !p.charsize
        charth = !p.charthick & xth = !x.thick & yth = !y.thick
        !p.thick = 2.0 & !p.charsize=1.5 & !p.charthick=2.0 & !x.thick=2 & !y.thick=2

        colortable1
        path = trgb_datapath()
        plotpath = path[3]+'SMF/' ; plot subdirectory

        if keyword_set(ps) then begin
            ps_open, plotpath+objname+'_smf', /ps_fonts
            device, /times, /inches
        endif else begin
            window, 0, xs=450, ys=450
            device, decomposed=0
        endelse

; ----------------------------------------------------------------------
; linear luminosity function
        plot, [min(lf.mag),max(lf.mag)], [0,max(phi)*1.1], /nodata, $
          line=0, yminor=2, xminor=3, $
          xsty=3, ysty=3, yrange = [0,max(phi)*1.1], $
          xtickname=replicate(' ',10), position=[0.12,0.32,0.504,0.92], color=1
        oplot, lf.mag, phi, color=3
        xyouts, [0.52,0.52], [0.94,0.94], infobase.truename, $
          /normal, charthick=2, charsize=2.2, align=0.5, color=1
        xyouts, [0.006,0.006], [0.62,0.62], 'N', $
          /normal, align=0.5, orient=90
; ----------------------------------------------------------------------
; linear TRGB
        oplot, [trgblin,trgblin], [!y.crange[0],!y.crange[1]], line=2, thick=3, color=2
; ----------------------------------------------------------------------
;; legend
;        legend, ['Stars '+strn(nstars), 'Binsize '+strn(lf.binsize,format='(F4.2)'), $
;                 'TRGB '+strn(trgblin,format='(F6.2)')+$
;                 ' '+strn(trgblinerr,format='(F4.2)')], $
;          charsize = 1.2, box=0 ;, /right, /bottom
; ----------------------------------------------------------------------
; linear response
        plot, [min(lf.mag),max(lf.mag)], $
          [min(resplin)>0.,max(resplin)*1.3], /nodata, $
          line=0, xsty=3, ysty=3, position=[0.12,0.12,0.504,0.32], /noerase, $
          xtit='I', yminor=2, xminor=3, color=1, yrange=[min(resplin)>0.,max(resplin)*1.3]
        oplot, lf.mag, resplin, color=5
        oplot, [trgblin,trgblin], [!y.crange[0],!y.crange[1]], line=2, thick=3, color=2
        xyouts, [0.006,0.006], [0.22,0.22], 'Linear Response', $
          /normal, align=0.5, orient=90, color=1
; ----------------------------------------------------------------------
; log luminosity function
        plot, [min(lf.mag),max(lf.mag)], [min(alog10(phi>1.)),max(alog10(phi>1.))*1.1], /nodata, $
          line=0, /noerase, yminor=2, xminor=3, xsty=3, ysty=7, $
          yrange = [min(alog10(phi>1.)),max(alog10(phi>1.))*1.1], position=[0.504,0.32,0.908,0.92], $
          ytickname = replicate(' ',10), xtickname=replicate(' ',10), color=1
        axis, yaxis=1, yminor=2, xminor=3, ysty=3, color=1
        oplot, lf.mag, alog10(phi>1.), color=7
        xyouts, [0.994,0.994], [0.62,0.62], 'log N', $
          /normal, align=0.5, orient=90
; ----------------------------------------------------------------------
; log TRGB
        oplot, [trgblog,trgblog], [!y.crange[0],!y.crange[1]], line=2, thick=3, color=2
; ----------------------------------------------------------------------
;; legend
;        legend, ['Stars '+strn(nstars), 'Binsize '+strn(lf.binsize,format='(F4.2)'), $
;                 'TRGB '+strn(trgblog,format='(F6.2)')+$
;                 ' '+strn(trgblogerr,format='(F4.2)')], $
;          charsize = 1.2, box=0 ;, /right, /bottom
; ----------------------------------------------------------------------
; log response
        plot, [min(lf.mag),max(lf.mag)], [min(resplog),max(resplog)*1.3], line=0, $
          xsty=3, ysty=7, position=[0.504,0.12,0.908,0.32], /noerase, /nodata, $
          ytickname=replicate(' ',10), yminor=2, xminor=3, xtit='I ', color=1, $
          yrange=[min(resplog),max(resplog)*1.3]
        axis, yaxis=1, yminor=2, xminor=3, ysty=3, color=1
        oplot, lf.mag, resplog, color=5
;       oploterror, lf.mag, resplog, respsiglog, errcolor=5
        xyouts, [0.994,0.994], [0.22,0.22], 'Log Response', $
          /normal, align=0.5, orient=90, color=1
        oplot, [trgblog,trgblog], [!y.crange[0],!y.crange[1]], line=2, thick=3, color=2
; ----------------------------------------------------------------------

        !x.ticklen=0
        !p.thick = p_thick & !p.charsize=chars & !p.charthick=charth
        !x.thick = xth & !y.thick=yth

        if keyword_set(ps) then ps_close
        
return
end
