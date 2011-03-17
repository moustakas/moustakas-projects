pro ages_check_speclinefit, _extra=extra
; jm04sep14uofa
    
    ages = read_ages(_extra=extra)
    ngalaxy = n_elements(ages)
    if (ngalaxy eq 0L) then return
    
    linename = strcompress(ages[0].linename,/remove)
    nline = n_elements(linename)

    tags = tag_names(ages)
    ntags = n_elements(tags)

    badindx = -1L
    
    for j = 0L, nline-1L do begin

       splog, 'Checking '+linename[j]+'.'

       keeptags = [linename[j]+['','_CONTINUUM','_EW','_SIGMA']]
       line = struct_trimtags(ages,select=keeptags)
       tags = tag_names(line)
;      line.fit_id = line.fit_id*3

       for ii = 0L, n_tags(line)-1L do begin
          inf = where((finite((line.(ii))[0,*]) eq 0B) or (finite((line.(ii))[1,*]) eq 0B),ninf)
          if (ninf eq 0L) then splog, 'No INFINITY/NAN detected for tag '+tags[ii] else begin
             splog, 'INFINITY/NAN.'
             struct_print, struct_trimtags(ages[inf],select=['GALAXY','Z'])
          endelse
       endfor
          
;      cc = get_kbrd(1)
       print

    endfor

return
end    
