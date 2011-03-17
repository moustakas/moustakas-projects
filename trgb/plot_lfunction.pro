pro plot_lfunction, name, lfunction1, lfunction2, log=log
;+
; NAME:
;	PLOT_LFUNCTION
;
; PURPOSE:
;	Plot a stellar luminosity function.
;
; INPUTS:
;	name	  : string name of the object
;	lfunction1 : luminosity function structure returned by trgb_lfunction
;
; OPTIONAL INPUTS:
;	lfunction2 : overplot a second luminosity function
;
; KEYWORD PARAMETERS:
;	log : plot a logarithmic luminosity function
;
; OUTPUTS:
;
; COMMON BLOCKS:
;
; MODIFICATION HISTORY:
;	John Moustakas, 2000 June 20, UCB
;	generalized, jm00aug4ucb
;-

        npar = n_params()
        if npar gt 2 then begin

            if keyword_set(log) then begin
                lf1 = alog10(float(lfunction1.hist)>1.)
                lf2 = alog10(float(lfunction2.hist)>1.)
                ytitle = 'log N(m)'
            endif else begin
                lf1 = lfunction1.hist
                lf2 = lfunction2.hist
                ytitle = 'N(m)'
            endelse
                
            ymax = max(lf1) > max(lf2)

        endif else begin

            if keyword_set(log) then begin
                lf1 = alog10(float(lfunction1.hist)>1.)
                ytitle = 'log N(m)'
            endif else begin
                lf1 = lfunction1.hist
                ytitle = 'N(m)'
            endelse

            ymax = max(lf1)

        endelse

        binsize = lfunction1.binsize
        xmag = lfunction1.mag
        nstars = n_elements(lfunction1.mags)
            
        plot, xmag, lf1, ps=10, line=0, yminor=3, color=7, $
          xsty=3, ysty=3, xminor=3, tit=name+' Luminosity Function', $
          ytitle=ytitle, xtit='I', yrange = [0,ymax]
        if npar gt 2 then oplot, xmag, lf2, ps=10, line=2, color=5 else begin
            legend, ['Stars '+strn(nstars), 'Binsize '+strn(lfunction1.binsize,format='(F5.2)')], $
              charsize = 1.5, textcolor=[1,1], box=0 ;, /right, /bottom
        endelse
            
return
end

