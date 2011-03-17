pro highzea_make_speclist, outfile, datapath=datapath
; jm06jun28uofa - based on ATLAS2D_MAKE_SPECLISTS

    if n_elements(outfile) eq 0L then begin
       print, 'Syntax - highzea_make_speclist, outfile, datapath='
       return
    endif
    
    if n_elements(datapath) eq 0L then datapath = cwd()
    pushd, datapath
    
    flist = file_search('caliblist*.txt',count=fcount)
    for i = 0L, fcount-1L do begin
       readcol, flist[i], filename, format='A', comment='#', /silent
       if i eq 0L then fitsfile = filename else fitsfile = [fitsfile,filename]
       
    endfor

    nfits = n_elements(fitsfile)
    object = strarr(nfits)
    for j = 0L, nfits-1L do object[j] = strupcase(sxpar(headfits(fitsfile[j]),'OBJECT'))

    len = string(max(strlen('fw'+fitsfile)))
       
    openw, lun1, outfile, /get_lun
    for k = 0L, nfits-1L do printf, lun1, 'fw'+fitsfile[k], object[k], object[k], $
      format='(A'+len+',A25,A25)'
    free_lun, lun1

    popd

return
end    
