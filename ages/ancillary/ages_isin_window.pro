function ages_isin_window, ra, dec
; jm10feb01ucsd - return a bit indicating whether or not an object is
;   in the AGES main survey area

    common ages_window, windowpoly
    
    nobj = n_elements(ra)
    if (nobj eq 0) or (n_elements(dec) eq 0) then begin
       doc_library, 'ages_isin_window'
       return, -1
    endif

    if (n_elements(windowpoly) eq 0) then begin
       polyfile = ages_path(/window)+'ages_window_final.ply'
       if (file_test(polyfile,/reg) eq 0) then begin
          splog, 'Polygon file '+polyfile+' not found!'
          return, 0.0
       endif
       read_mangle_polygons, polyfile, windowpoly
    endif

    yesno = is_in_window(ra=ra,dec=dec,windowpoly)

return, yesno
end
