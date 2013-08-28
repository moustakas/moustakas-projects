pro mzoplot_kk04_models, model_logr23=model_logr23, model_logq=model_logq, $
  linestyle=linestyle, linecolor=linecolor, nolegend=nolegend, _extra=extra
    if (n_elements(linestyle) eq 0) then linestyle = [0,5,3,4,2]
    if (n_elements(linecolor) eq 0) then linecolor = $
      replicate('default',n_elements(linestyle))
    if (n_elements(model_logq) eq 0) then $
      model_logq = alog10([1.2D8,8D7,4D7,2D7,1D7]) ; alog10(4E7)
    if (n_elements(model_logr23) eq 0) then model_logr23 = range(-0.5,1.1,1500)
    for iq = 0L, n_elements(model_logq)-1L do begin
       model_logoh_upper = 9.72D - 0.777*model_logr23 - $
         0.951*model_logr23^2 - 0.072*model_logr23^3 - $
         0.811*model_logr23^4 - model_logq[iq]*(0.0737 - $
         0.0713*model_logr23 - 0.141*model_logr23^2 + $
         0.0373*model_logr23^3 - 0.058*model_logr23^4)
       model_logoh_lower = 9.40D + 4.65D*model_logr23 - $
         3.17D*model_logr23^2 - model_logq[iq]*$
         (0.272D + 0.547D*model_logr23 - 0.513D*model_logr23^2)
       model_good1 = where((model_logoh_upper gt model_logoh_lower))
       model_good2 = where((model_logoh_lower[model_good1] gt 7.5))
       djs_oplot, model_logr23[model_good1], model_logoh_upper[model_good1], $
         linestyle=linestyle[iq], color=linecolor[iq], _extra=extra
       djs_oplot, model_logr23[model_good1], model_logoh_lower[model_good1], $
         linestyle=linestyle[iq], color=linecolor[iq], _extra=extra
    endfor
    if (keyword_set(nolegend) eq 0) then begin
       legend, 'log(U)='+string(model_logq-alog10(im_light(/cm)),format='(F5.2)'), $
         /right, /bottom, box=0, charsize=1.3, line=linestyle, pspacing=1.4
    endif
return
end
