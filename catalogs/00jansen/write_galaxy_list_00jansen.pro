pro write_galaxy_list_00jansen, data, write=write
; jm04nov29uofa - generate a list of galaxy names that are NED-compatible

    path = getenv('CATALOGS_DIR')+'/00jansen/'

    data = rsex(path+'nfgs_sample_properties.dat')
    ngalaxy = n_elements(data)

    this = where((strmatch(data.name,'*NGC5491*') eq 1B),nthis)
    if (nthis ne 0L) then data[this].ugc = 'NGC5491'
    
    mrk = where((strmatch(data.name,'*A12195+7535*') eq 1B),nmrk)
    if (nmrk ne 0L) then data[mrk].ugc = 'MRK0205'
    
    mrk2 = where((strmatch(data.name,'*A15016+1037*') eq 1B),nmrk2)
    if (nmrk2 ne 0L) then data[mrk2].ugc = 'MRK0841'

    if keyword_set(write) then begin
       openw, lun, path+'galaxy_list.txt', /get_lun
       niceprintf, lun, data.ugc
       free_lun, lun
    endif    

return
end    
