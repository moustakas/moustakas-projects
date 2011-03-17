pro plot_gdds, file, _extra=extra
; jm04mar28uofa
; plot data from the Gemini Deep Deep Survey

    readcol, file, wave, flux, ferr, format='F,F,F', comment='#', $
      /silent, delimiter=' '

    djs_plot, wave, flux, ps=10, xrange=[6500,9800], xsty=3, ysty=3, $
      xthick=2.0, ythick=2.0, charsize=1.5, charthick=2.0

return
end
