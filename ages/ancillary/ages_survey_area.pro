function ages_survey_area
; jm10feb01ucsd - read the polygons written out by AGES_WINDOW to get
; the total non-overlapping survey area for AGES

    polyfile = ages_path(/window)+'ages_window_final.ply'
    if (file_test(polyfile,/reg) eq 0) then begin
       splog, 'Polygon file '+polyfile+' not found!'
       return, 0.0
    endif

    read_mangle_polygons, polyfile, windowpoly

    area = 0.0D
    for ii = 0, n_elements(windowpoly)-1 do $
      area = area + garea(windowpoly[ii])

return, area*(180.0D/!dpi)^2.0 ; [deg^2]
end
