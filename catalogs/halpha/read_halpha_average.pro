function read_halpha_average
; jm04may04uofa
; read Kennicutt's average literature H-alpha catalog

    hapath = filepath('',root_dir=getenv('CATALOGS_DIR'),subdirectory='halpha')

    pushd, hapath

    catfile = 'halpha_average.dat'
    if file_test(catfile,/regular) eq 0L then begin

       splog, 'Catalog '+hapath+catfile+' not found.'
       return, -1L

    endif

    data = djs_readlines(catfile,nhead=7)
    ngalaxy = n_elements(data)

    cat = {$
      name:          ' ', $
      ra:            ' ', $
      dec:           ' ', $
      ew:          -999D0,$
      ha_flux:     -999D0}
    cat = replicate(cat,ngalaxy)

    for k = 0L, ngalaxy-1L do begin

       line = data[k]
       
       cat[k].name = strcompress(repstr(strmid(line,0,17),'.',''),/remove)
       cat[k].ra = repstr(repstr(repstr(strmid(line,17,13),'h',':'),'m',':'),'s',':')
       cat[k].dec = repstr(repstr(repstr(strmid(line,30,12),'h',':'),'m',':'),'s',':')

       ew = double(strmid(line,42,9))
       flux = double(strmid(line,51,7))

       if (ew ne 0.0) and (ew gt -90.0) then cat[k].ew = ew
       if (flux ne 0.0) then cat[k].ha_flux = -flux
       
    endfor

;   struct_print, cat

    popd

return, cat
end
