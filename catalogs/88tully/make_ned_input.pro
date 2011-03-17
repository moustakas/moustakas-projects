pro make_ned_input

    path = '/home/ioannis/idl/catalogs/88tully'
    pushd, path
    tully = read_88tully()
    ntully = n_elements(tully)

    search = 5.0

    ra = tully.x_raj2000
    de = tully.x_dej2000
    
    openw, lun, 'ned_input.txt', /get_lun
    for i = 0L, ntully-1L do printf, lun, strmid(ra[i],0,2)+'h'+strmid(ra[i],3,2)+'m ; '+$
      strmid(de[i],0,3)+'d'+strmid(de[i],4,2)+'m ; ', strn(search), format='(A20,A8)'
    free_lun, lun

    popd

return
end    
