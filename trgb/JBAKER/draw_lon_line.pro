
PRO draw_lon_line, $
                   lon, $
                   NPOINTS=npts, $
                   THICK=thick, $
                   LINESTYLE=linestyle

IF NOT keyword_set( npts ) THEN npts=500
if NOT keyword_set( thick ) THEN thick=0.5*!p.thick

lat = 180. * findgen(npts)/float(npts-1) - 90.
plots, fltarr(npts)-lon, lat, thick=thick, linestyle=linestyle

end

