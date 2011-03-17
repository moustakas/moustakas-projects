function read_halpha_literature
; jm03may24uofa
; read Kennicutt's literature H-alpha catalog

    hapath = filepath('',root_dir=getenv('CATALOGS_DIR'),subdirectory='halpha')

    pushd, hapath

    catfile = 'halpha.dat'
    if file_test(catfile,/regular) eq 0L then begin

       splog, 'Catalog '+hapath+catfile+' not found.'
       return, -1L

    endif

    data = djs_readlines(catfile,nhead=1)
    ngalaxy = n_elements(data)

    cat = {$
      name1:         ' ', $
      name2:         ' ', $
      ra:            ' ', $
      dec:           ' ', $
      A_B:         -999D0,$
      type:          ' ', $
      B:           -999D0,$
      dmaj:        -999D0,$
      dmin:        -999D0,$
      cz:          -999D0,$
      aperture:       '', $
      ew:          -999D0,$
      ew_err:      -999D0,$
      ha_flux:     -999D0}
    cat = replicate(cat,ngalaxy)

    for k = 0L, ngalaxy-1L do begin

       line = data[k]
;      len = strlen(line)
;      for j = 0L, len-1L do niceprint, j, strmid(line,j,1)
       
       cat[k].name1       = strcompress(strmid(line,0,19),/remove)
       cat[k].name2       = strcompress(strmid(line,19+1,22),/remove)
       cat[k].ra          = strjoin(strsplit(im_dec2hms(im_hms2dec(strmid(line,19+22+1,15))),' ',/extract),':')
       cat[k].dec         = strjoin(strsplit(im_dec2hms(im_hms2dec(strmid(line,19+22+15+1,14))),' ',/extract),':')
       cat[k].A_B         = strcompress(strmid(line,19+22+15+14+1,4),/remove)
       cat[k].type        = strmid(line,19+22+15+14+4+1,18)

       B = strmid(line,19+22+15+14+4+18+1,6)
       if strcompress(B,/remove) ne '' then cat[k].B = B

       dmaj = strmid(line,19+22+15+14+4+18+6+1,5)
       if strcompress(dmaj,/remove) ne '' then cat[k].dmaj = dmaj
       
       dmin = strmid(line,19+22+15+14+4+18+6+5+1,5)
       if strcompress(dmin,/remove) ne '' then cat[k].dmin = dmin

       cz = strmid(line,19+22+15+14+4+18+6+5+5+1,7)
       if strcompress(cz,/remove) ne '' then cat[k].cz = cz
       
       space = 19
       
       ew = strmid(line,19+22+15+14+4+18+6+5+5+7+space+1,6)
       if strcompress(ew,/remove) ne '' then cat[k].ew = ew

       ew_err = strmid(line,19+22+15+14+4+18+6+5+5+7+space+6+1,4)
       if strcompress(ew_err,/remove) ne '' then cat[k].ew_err = ew_err

       ha_flux = strmid(line,19+22+15+14+4+18+6+5+5+7+space+6+4+1,6)
       if (strcompress(ha_flux,/remove) ne '') then begin
          if (ha_flux ne 0.0) then cat[k].ha_flux = -ha_flux
       endif

;      if strmatch(cat[k].name1,'*4450*') then stop
       
    endfor

;   struct_print, cat

    popd

return, cat
end
