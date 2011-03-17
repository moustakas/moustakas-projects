pro plot_response, name, resplin, respsiglin, resplog, respsiglog, xmag
;+
; NAME:
;	PLOT_RESPONSE
;
; PURPOSE:
;	Plot the Sobol response filter.
;
; INPUTS:
;	name 	   : galaxy name
;	resplin    : Sobol response of a linear luminosity function
;	respsiglin : uncertainty in resplin
;	resplog	   : Sobol response of a logarithmic luminosity
;	     	     function
;	respsiglog : uncertainty in resplog
;	xmag	   : luminosity function ordinate (magnitudes)
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;
; PROCEDURES USED:
;	
; MODIFICATION HISTORY:
;	John Moustakas, 2000 June 26, UCB
;-

; trap errors

	goodlin = where(finite(respsiglin) eq 1L)
	goodlog = where(finite(respsiglog) eq 1L)

; restrict the magnitude range for plotting purposes

        xlinminmax = minmax(xmag[goodlin])
        xlogminmax = minmax(xmag[goodlog])
        
        !p.multi = [0,1,2]
        plot, xmag[goodlin], resplin[goodlin], line=0, color=5, ps=10, $
          xsty=1, ysty=1, tit=name+' Linear Response', $
          ytit = 'Response', xtit='I Magnitude', yrange=[0,max(resplin[goodlin])*1.2], $
          xrange = [xlinminmax[0],xlinminmax[1]]
        oploterror, xmag[goodlin], resplin[goodlin], respsiglin[goodlin], errcolor=5
        oplot, [!x.crange[0],!x.crange[1]], [0,0], color=1, thick=1.5
        plot, xmag, resplog, line=0, color=3, ps=10, $
          xsty=1, ysty=1, tit=name+' Log Response', $
          ytit = 'Response', xtit='Magnitude', yrange=[0,max(resplog[goodlog])*1.2], $
          xrange = [xlogminmax[0],xlogminmax[1]]
        oploterror, xmag[goodlog], resplog[goodlog], respsiglog[goodlog], errcolor=3
        oplot, [!x.crange[0],!x.crange[1]], [0,0], color=1, thick=1.5
        !p.multi = 0

return
end

