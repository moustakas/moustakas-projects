pro oplot_lzfit, coeff, band=band, linestyle=linestyle, linecolor=linecolor, $
  thick=thick, tremonti=tremonti, evol=evol, _extra=extra

    pivotmag = lz_pivotmag()
    
    ohrange = !y.crange
    ohoffset = 0.07 ; 1
    magaxis = range(!x.crange[0]-0.3,!x.crange[1]+0.3,50)

    if keyword_set(tremonti) then begin
       lz = poly(magaxis,[5.283,-0.185])
       inrange = where((lz gt ohrange[0]+ohoffset) and $
         (lz lt ohrange[1]-ohoffset))
       djs_oplot, magaxis[inrange], lz[inrange], thick=10, $
         line=linestyle, color=im_color(linecolor)
       return
    endif

    if keyword_set(evol) then $
      ohmodel = lzevol_func(magaxis,coeff,pivotmag=pivotmag,_extra=extra) else $
        ohmodel = poly(magaxis-pivotmag,coeff)
    
    inrange = where((ohmodel gt ohrange[0]+ohoffset) and $
      (ohmodel lt ohrange[1]-ohoffset))
    djs_oplot, magaxis[inrange], ohmodel[inrange], $
      thick=7, linestyle=linestyle, color=im_color(linecolor)

return
end
