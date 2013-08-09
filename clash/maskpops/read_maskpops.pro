function read_maskpops1, file

    nfile = n_elements(file)
    if (nfile eq 0) then begin
       splog, 'Need FILE input'
       return, -1
    endif

; call recursively
    if (nfile gt 1) then begin
       for ii = 0, nfile-1 do begin
          cat1 = read_maskpops1(file[ii])
          if (ii eq 0) then cat = cat1 else cat = [cat,cat1]
       endfor
    endif else begin
       if file_test(file) eq 0 then message, 'No file found!'
       phot = rsex(file)
;      niceprint, phot.abmag_zpt & print
       filt = maskpops_filterlist(short=short)
;      filt = clash_filterlist(short=short)
       nfilt = n_elements(filt)

       cat = {prefix: '', z: 0.0}
       for ii = 0, nfilt-1 do cat = create_struct(cat,$
         short[ii]+'_flux',float(phot[ii].fnu),short[ii]+'_fluxerr',$
         float(phot[ii].fnuerr))
    endelse

return, cat
end

function read_maskpops
; jm12may04ucsd - read the MASKPOPS output catometric catalog and put
; it into the standard iSEDfit catalog format

    path = maskpops_path()
;   metafile = path+'roiphot_arcs.dat'
    metafile = path+'roiphot_startup.dat'
    meta = rsex(metafile)

;   file = strtrim(meta.rootname,2)+'_'+strtrim(meta.scaledir,2)+'_phot.cat'
    file = strtrim(meta.cluster,2)+'_'+strtrim(meta.rootname,2)+'_'+$
      strtrim(meta.scaledir,2)+'_'+strtrim(meta.datestamp,2)+'_phot.cat'
    cat = read_maskpops1(path+file)

    cat.z = meta.redshift
    cat.prefix = strtrim(meta.rootname,2)

; we *need* to sort by redshift here, otherwise iSEDfit's
; USE_REDSHIFT chokes
    cat = cat[sort(cat.z)]
    
return, cat
end
