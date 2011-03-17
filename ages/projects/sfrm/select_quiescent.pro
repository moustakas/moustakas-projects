function select_quiescent, nuvmr, rmj, active=active, pars=pars
; jm10feb06ucsd - select quiescent and actively star-forming galaxies
;   based on the NUV-r vs r-J color-color diagram (see Williams et
;   al. 2009, Ilbert et al. 2010); NUVMR = NUV-r, rmj = r-J

    pars = {y1: 4.2}
    if (n_elements(nuvmr) gt 0) then begin
       quiescent = where((nuvmr gt pars.y1),nquiescent,comp=active)
    endif else quiescent = -1

;; ## old three-sided selection box a la Williams+09
;    pars = {$
;      x1: 0.0, x2: 0.9, x3: 1.7, $
;      y1: 3.5, y2: 5.2, y3: 7.0}
;
;    slope = (pars.y2-pars.y1)/(pars.x3-pars.x2)
;    int = pars.y1-slope*pars.x2
;
;    if (n_elements(nuvmr) gt 0) then begin
;       quiescent = where(((rmj lt pars.x2) and (nuvmr gt pars.y1)) or $
;         ((rmj gt pars.x2) and (rmj lt pars.x3) and $
;         (nuvmr gt poly(rmj,[int,slope]))),nquiescent,comp=active)
;    endif else quiescent = -1

return, quiescent
end
