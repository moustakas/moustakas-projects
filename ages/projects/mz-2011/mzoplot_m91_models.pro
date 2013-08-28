pro mzoplot_m91_models, model_logr23=model_logr23, model_logo32=model_logo32, $
  linestyle=linestyle, linecolor=linecolor, nolegend=nolegend, _extra=extra
    if (n_elements(linestyle) eq 0) then linestyle = [0,5,3,4,2]
    if (n_elements(linecolor) eq 0) then linecolor = $
      replicate('default',n_elements(linestyle))
    if (n_elements(model_logo32) eq 0) then $
      model_logo32 = range(-1.0,0.2,5)
    if (n_elements(model_logr23) eq 0) then model_logr23 = range(-0.5,1.1,1500)
    for ii = 0L, n_elements(model_logo32)-1L do begin
       model_logoh_lower = (12.0 - 4.944 + 0.767*model_logr23 + 0.602*model_logr23^2) - $
         model_logo32[ii]*(0.29 + 0.332*model_logr23 - 0.331*model_logr23^2)
       model_logoh_upper = (12.0 - 2.939 - 0.2*model_logr23 - 0.237*model_logr23^2 - $
         0.305*model_logr23^3 - 0.0283*model_logr23^4) - model_logo32[ii]*(0.0047 - $
         0.0221*model_logr23 - 0.102*model_logr23^2 - 0.0817*model_logr23^3 - 0.00717*model_logr23^4)
       model_good1 = where((model_logoh_upper gt model_logoh_lower))
       model_good2 = where((model_logoh_lower[model_good1] gt 7.5))
       djs_oplot, model_logr23[model_good1], model_logoh_upper[model_good1], $
         linestyle=linestyle[ii], color=linecolor[ii], _extra=extra
       djs_oplot, model_logr23[model_good1], model_logoh_lower[model_good1], $
         linestyle=linestyle[ii], color=linecolor[ii], _extra=extra
    endfor
    if (keyword_set(nolegend) eq 0) then begin
       legend, textoidl('log(O_{23})='+string(model_logo32,format='(F5.2)')), $
         /right, /bottom, box=0, charsize=1.3, line=linestyle, pspacing=1.4
    endif
return
end
