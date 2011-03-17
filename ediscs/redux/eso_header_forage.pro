function eso_header_forage, h
; jm04jun20uofa
; read an ESO non-standard header

    good = where(strcompress(h,/remove) ne '',ngood)
    if (ngood ne 0L) then hmine = h[good] else hmine = h

    nkeys = n_elements(hmine)
    for i = 0L, nkeys-1L do begin

       info = strtrim(strsplit(hmine[i],'=/',/extract),2)

       if (strmatch(hmine[i],'*comment*',/fold) eq 0B) and $
         (strmatch(hmine[i],'*history*',/fold) eq 0B) and $
         (strmatch(hmine[i],'*end*',/fold) eq 0B) then begin
       
          tag = repstr(repstr(repstr(info[0],'-','_'),' ','_'),'HIERARCH_ESO_','')
          value = repstr(info[1],"'",'')
       
          if (n_elements(f) eq 0L) then f = create_struct(tag,value) else $
            f = create_struct(f,create_struct(tag,value))

       endif

    endfor

return, f
end
