pro build_deep2_spec1d_tarball
; jm13jul15siena - make a tarball of all the 1D spectra in
; SPECLINEFIT.GOODSPEC1D.Q34.FITS

    path = deep2_path(/dr4)
;   spawn, 'mkdir -f '+path+'tmp/'
    
    jj = read_deep2(/ispec)
;   for ii = 0, n_elements(jj)-1 do file_copy, path+strtrim(jj[ii].file,2), path+'/tmp/';, /overwrite
;   for ii = 0, 3 do file_copy, path+strtrim(jj[ii].file,2), '/tmp/', /overwrite

    pushd, path
    files = strtrim(jj.file,2)
;   files = strtrim(file_basename(jj.file),2)
    spawn, 'tar -cvf '+path+'deep2.dr4.spec1d.tar '+files[0], /sh
;   for ii = 1, 10 do spawn, 'tar -rvf '+path+'deep2.dr4.spec1d.tar '+files[ii], /sh
    for ii = 1, n_elements(jj)-1 do spawn, 'tar rvf '+path+'deep2.dr4.spec1d.tar '+files[ii], /sh
    spawn, 'gzip -f '+path+'deep2.dr4.spec1d.tar', /sh
    popd

stop    
    
    
return
end
    
