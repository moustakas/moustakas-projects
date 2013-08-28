pro mzplot_scatterplot, x, y, sdss=sdss, nooutliers=nooutliers, $
  levels=levels, annotate=annotate, npix=npix, outpsym=outpsym, $
  outsymsize=outsymsize, outcolor=outcolor, _extra=extra
; mass-metallicity plot wrapper on HOGG_SCATTERPLOT
    if keyword_set(sdss) then begin
       if (n_elements(levels) eq 0) then $
         levels = [0.01,0.05,0.25,0.5,0.75,0.95,0.99]
       if (n_elements(npix) eq 0) then npix = 50
       if (n_elements(outsymsize) eq 0) then outsymsize = 0.15
    endif else begin
       if (n_elements(levels) eq 0) then $
         levels = [0.1,0.25,0.5,0.75,0.95]
       if (n_elements(npix) eq 0) then npix = 16
       if (n_elements(outsymsize) eq 0) then outsymsize = 0.4
    endelse
    if keyword_set(annotate) then cannotation = $
      string(levels,format='(F4.2)')
    if (n_elements(outpsym) eq 0) then outpsym = 16
    if (n_elements(outcolor) eq 0) then outcolor = im_color('grey60',101)
    hogg_scatterplot, x, y, outliers=(keyword_set(nooutliers) eq 0), $
      outpsym=symcat(outpsym), outsymsize=outsymsize, $
      outcolor=outcolor, _extra=extra, exp=0.5, xnpix=npix, $
      ynpix=npix, levels=levels, cannotation=cannotation, /internal
return
end
