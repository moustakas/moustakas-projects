pro update_ibatch, ibatch=ibatch, write=write
; jm08jan17nyu - NO LONGER NEEDED
    
    allpath = atlas_path(/spec2datlas)+[$
;     '94nov/',$ ; Turner
;     '95mar/',$ ; Turner
;     '95oct/',$ ; Turner
;     '96apr/',$ ; Turner
;     '97apr/',$ ; Turner
      '98mar/',$ ; done
;     '98apr/',$ ; do it by hand
      '98jun/',$
      '98oct/',$
      '99apr/',$
      '99may/',$
      '99nov/',$
;     '00apr/',$ ; do it by hand
;     '01nov/',$ ; done by hand
;     '01dec/',$ ; special case
      '02feb/',$
      '02apr/',$
      '02may/',$
      '03may/']
;     '05apr/',$ ; already done
;     '06mar/',$ ; done
;     '06may/']  ; done

;   path = atlas_path(/spec2datlas)+'00apr/'

    for ip = 0L, n_elements(allpath)-1L do begin

       path = allpath[ip]
       
       if keyword_set(ibatch) then begin
          batch = file_search(path+'ibatch_*.txt',count=cc)
          for ii = 0L, cc-1L do begin
; add TRACENAME to the ibatch files
             suffix = strmid(file_basename(batch[ii]),7,7)
             file = djs_readlines(batch[ii])
             file = repstr(file,'TRACENAME                               ',$
               'TRACENAME          '+'trace_'+suffix+'.idlsave')
             if keyword_set(write) then begin
                splog, 'Re-writing '+batch[ii]
                openw, lun, batch[ii], /get_lun
                niceprintf, lun, file
                free_lun, lun
             endif else begin
                print, file_basename(batch[ii])+': '+$
                  file[where(strmatch(file,'*TRACENAME*'))]
;            hgrep, file, 'TRACENAME'
;            niceprint, file
             endelse
          endfor
; run ISPEC_MAKELISTS
          flist = file_search(path+'a.*.fits')
          flist = strmid(file_basename(flist),0,4)
          root = flist[uniq(flist,sort(flist))]
          root = root[where((strmatch(root,'a.00*') eq 0B))]

          prolist = file_search(path+'*_script.pro')
          spawn, 'grep "dateroot = " '+prolist, daterootline
          spawn, 'grep "dates = " '+prolist, datesline
          junk = execute(daterootline[0])
          junk = execute(datesline[0])
          if keyword_set(write) then begin
             splog, 'Running ispec_makelists'
             pushd, path
             call_procedure, 'ispec_makelists', root, dates, /trace, /overwrite
             popd
          endif else begin
;            print, 'ROOT, DATES'
             niceprint, root, dates
          endelse
          if keyword_set(write) then begin
             print
;            cc = get_kbrd(1)
          endif else begin
             print, '-------------------------'
             cc = get_kbrd(1)
          endelse
       endif 
    endfor
       
return
end
    
