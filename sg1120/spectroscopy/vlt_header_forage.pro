function vlt_header_forage, h, silent=silent
; jm04nov05uofa
; read a VLT non-standard header

    good = where(strcompress(h,/remove) ne '',ngood)
    if (ngood ne 0L) then hmine = h[good] else hmine = h

    nkeys = n_elements(hmine)
    for i = 0L, nkeys-1L do begin

       info = strtrim(strsplit(hmine[i],'=/',/extract),2)

       if (strmatch(hmine[i],'*comment*',/fold) eq 0B) and $
         (strmatch(hmine[i],'*history*',/fold) eq 0B) and $
         (strmatch(info[0],'END') eq 0B) then begin
       
          tag = repstr(repstr(repstr(info[0],'-','_'),' ','_'),'HIERARCH_ESO_','')
          value = repstr(info[1],"'",'')

          if (n_elements(f) eq 0L) then begin

             f = create_struct(tag,value)

          endif else begin

             exist = where(strmatch(tag_names(f),tag,/fold) eq 1B,nexist)
             if (nexist ne 0L) then begin
                if not keyword_set(silent) then splog, 'Duplicate structure tag '+tag+'.' 
             endif else f = create_struct(f,create_struct(tag,value))

          endelse 

       endif 

    endfor

return, f
end
