pro oplot_mzfit, coeff, pivotmass, bin_mass=bin_mass, bin_oh=bin_oh, $
  linestyle=linestyle, linecolor=linecolor, $
  thick=thick, psymcolor=psymcolor, nopsym=nopsym, tremonti_plot=tremonti_plot, $
  tremonti_linecolor=tremonti_linecolor, tremonti_linestyle=tremonti_linestyle, $
  tremonti_extrap=tremonti_extrap, massaxis=massaxis, mshift=mshift

; overplot the fitted MZ relation
    if keyword_set(tremonti_plot) then begin    
       tmassaxis = im_array(8.8,11.3,0.01) ; Kroupa+01 IMF
       mztremonti = tremonti_mz(massaxis=tmassaxis,pivotmass=pivotmass,/kk04)
       djs_oplot, tmassaxis, mztremonti, line=tremonti_linestyle, $
         thick=10.0, color=fsc_color(tremonti_linecolor,99)
    endif
    if (keyword_set(nopsym) eq 0) then begin
       ww = where(bin_mass gt 0.0)
       djs_oplot, bin_mass[ww], bin_oh[ww], psym=symcat(6,thick=7), $
         symsize=0.6, color=fsc_color(psymcolor,100)
    endif
    if (n_elements(massaxis) eq 0L) then $
      massaxis = im_array(8.7,11.3,0.01)
    if (n_elements(mshift) eq 0L) then mshift = 0.0
    djs_oplot, massaxis+mshift, poly(massaxis-pivotmass,coeff), $
      thick=10.0, linestyle=linestyle, color=fsc_color(linecolor,101)

return
end
    
