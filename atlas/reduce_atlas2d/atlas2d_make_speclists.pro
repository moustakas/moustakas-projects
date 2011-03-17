pro atlas2d_make_speclists, outfile, datapath=datapath
; jm02oct9uofa
; generate a template SPECLIST.TXT or ATLASLIST.TXT or SINGSLIST.TXT
; for output to the "final product" directory (with MAKE_ATLAS2D)
; based on the files in caliblist_*.txt.

    if n_elements(outfile) eq 0L then begin
       print, 'Syntax - make_speclist, outfile, datapath='
       return
    endif
    
    if n_elements(datapath) eq 0L then datapath = cwd()
    pushd, datapath
    
    flist = findfile('caliblist*.txt',count=fcount)
    for i = 0L, fcount-1L do begin
       readcol, flist[i], filename, format='A', /silent
       if i eq 0L then fitsfile = filename else fitsfile = [fitsfile,filename]
       
    endfor

    nfits = n_elements(fitsfile)
    object = strarr(nfits)
    for j = 0L, nfits-1L do object[j] = strupcase(sxpar(headfits(fitsfile[j]),'OBJECT'))

    len = string(max(strlen('fdw'+fitsfile)))
       
    openw, lun1, outfile, /get_lun
    for k = 0L, nfits-1L do printf, lun1, 'fdw'+fitsfile[k], object[k], object[k], $
      format='(A'+len+',A25,A25)'
    free_lun, lun1

    popd

return
end    
