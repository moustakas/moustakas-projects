pro oplot_select_quiescent
; jm10feb06ucsd - overplot the box for selecting quiescent galaxies

    junk = select_quiescent(pars=pars)
;   color = 'blue'
    line = 5
    thick = 9

    djs_oplot, !x.crange, [pars.y1,pars.y1], $
      line=line, color=color, thick=thick

;; ## old three-sided selection box a la Williams+09
;   djs_oplot, [pars.x1,pars.x2], [pars.y1,pars.y1], $
;     line=line, color=color, thick=thick
;   djs_oplot, [pars.x2,pars.x3], [pars.y1,pars.y2], $
;     line=line, color=color, thick=thick
;   djs_oplot, [pars.x3,pars.x3], [pars.y2,pars.y3], $
;     line=line, color=color, thick=thick

return
end
    
