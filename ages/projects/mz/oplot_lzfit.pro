pro oplot_lzfit, coeff, band=band, linestyle=linestyle, linecolor=linecolor, $
  thick=thick, tremonti_plot=tremonti_plot, tremonti_linecolor=tremonti_linecolor, $
  tremonti_linestyle=tremonti_linestyle, evol=evol, _extra=extra

    pivotmag = lz_pivotmag()
    
    ohrange = !y.crange
    ohoffset = 0.07 ; 1

    case band of
       'B': magaxis = im_array(-24.0,-15.0,0.05)
       'g': magaxis = im_array(-24.0,-15.0,0.05)
       'r': magaxis = im_array(-24.0,-15.0,0.05)
    endcase
    magaxis = reverse(magaxis)

    if keyword_set(tremonti_plot) then begin
       glztremonti = tremonti_glz(magaxis=magaxis,kk04=1)
       inrange = where((glztremonti gt ohrange[0]+ohoffset) and $
         (glztremonti lt ohrange[1]-ohoffset))
       djs_oplot, magaxis[inrange], glztremonti[inrange], $
         line=tremonti_linestyle, thick=10.0, color=fsc_color(tremonti_linecolor,99)
    endif

    if keyword_set(evol) then $
      ohmodel = lzevol_func(magaxis,coeff,pivotmag=pivotmag,_extra=extra) else $
        ohmodel = poly(magaxis-pivotmag,coeff)
    
    inrange = where((ohmodel gt ohrange[0]+ohoffset) and $
      (ohmodel lt ohrange[1]-ohoffset))
    djs_oplot, magaxis[inrange], ohmodel[inrange], $
      thick=10.0, linestyle=linestyle, color=fsc_color(linecolor,101)

return
end
