pro stellarlocus_decals_scatterplot, xall, yall, xx, yy, xk, yk, $
  xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange, $
  type=type, feh=feh, label=label, position=position

; I,II = supergiant; III,IV = giant; V = dwarf    
;   dwarf = where(strmatch(strtrim(type,2),'*V*') or $
;     strmatch(strtrim(type,2),'*IV*'))
;   giant = where(strmatch(strtrim(type,2),'*I*') and $
;     (strmatch(strtrim(type,2),'*IV*') eq 0))

    poordwarf = where(strmatch(strtrim(type,2),'*V*') and $
      (strmatch(strtrim(type,2),'*IV*') eq 0) and (feh lt 0.0))
    richdwarf = where(strmatch(strtrim(type,2),'*V*') and $
      (strmatch(strtrim(type,2),'*IV*') eq 0) and (feh ge 0.0))
    poorgiant = where(strmatch(strtrim(type,2),'*I*') and (feh lt 0.0))
    richgiant = where(strmatch(strtrim(type,2),'*I*') and (feh ge 0.0))

    hogg_scatterplot, xall, yall, xtitle=xtitle, ytitle=ytitle, $
      xrange=xrange, yrange=yrange, /xsty, /ysty, xthick=3.0, ythick=3.0, $
      charsize=2.0, charthick=3.0, /outliers, outcolor=djs_icolor(''), $
      /internal, position=position
    im_legend, ['Metal-Rich Giant','Metal-Poor Giant','Metal-Rich Dwarf',$
      'Metal-Poor Dwarf','Kurucz Models'], /right, /bottom, box=0, $
      charsize=1.6, charthick=3.0, psym=[1,1,6,6,4], thick=5.0, $
      color=djs_icolor(['red','blue','dark green','cyan blue',''])
    im_legend, textoidl(label), /left, /top, box=0, charsize=2.0, charthick=3.0
;   legend, ['Pickles Giant','Pickles Dwarf','Kurucz Models'], /right, $
;     /bottom, box=0, charsize=1.5, charthick=3.0, psym=[1,6,4], thick=5.0, $
;     color=djs_icolor(['red','dark green','blue'])
    djs_oplot, xk, yk, psym=4, color='', thick=5.0
;   djs_oplot, xx[giant], yy[giant], psym=1, color='red', thick=5.0
;   djs_oplot, xx[dwarf], yy[dwarf], psym=6, color='dark green', thick=5.0
    djs_oplot, xx[richgiant], yy[richgiant], psym=1, symsize=1.3, color='red', thick=5.0
    djs_oplot, xx[poorgiant], yy[poorgiant], psym=1, symsize=1.3, color='blue', thick=5.0
    djs_oplot, xx[richdwarf], yy[richdwarf], psym=6, symsize=1.0, color='dark green', thick=5.0
    djs_oplot, xx[poordwarf], yy[poordwarf], psym=6, symsize=1.0, color='cyan blue', thick=5.0
    
return
end

