; jm00oct19
; remove duplicate points in the sed files
pro fix_seds

    on_error, 2

    path = '/home/ioannis/sirtf/sed/'
    flist = findfile(path)
    good = where(strpos(flist,'.txt') gt -1L,fcount)
    sedfiles = flist[good]

    openw, lun, path+'fix_seds.log', /get_lun ; write a log file of the changes
    printf, lun, '#'
    printf, lun, '# Entries at the following wavelengths were removed from the following SEDs:'
    printf, lun, '#'

    for j = 0L, fcount-1L do begin ; loop on the sed files

       spawn, ['\cp '+path+sedfiles[j]+' '+path+'OLD/'+sedfiles[j]]
       readcol, path+sedfiles[j], lambda, lum, /silent, format='D,D'

       nlambda = n_elements(lambda)
       nlum = n_elements(lum)

       if nlambda ne nlum then message, "There's something wrong with the SED."

       printf, lun, sedfiles[j]
       printf, lun, ' '

       openw, lun2, path+sedfiles[j], /get_lun ; overwrite the sed
       printf, lun2, lambda[0], lum[0], format='(2x,E9.3,2x,E10.4)'
       for k = 1L, nlambda-1L do begin
          if lambda[k] ne lambda[k-1L] then $
            printf, lun2, lambda[k], lum[k], format='(2x,E9.3,2x,E10.4)' else $
            printf, lun, lambda[k], lum[k], format='(2x,E9.3,2x,E10.4)'
       endfor
       printf, lun, ' '
       free_lun, lun2

    endfor
    free_lun, lun

return
end
