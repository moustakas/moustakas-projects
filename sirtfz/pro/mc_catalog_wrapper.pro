pro mc_catalog_wrapper
; jm01sep29uofa

    f1 = ['IRAC 8','MIPS 24','MIPS 70','MIPS 160']
    f2 = ['R DWFS','IRAC 3.6','MIPS 70']
        
    filtarray = ptrarr(2)
    filtarray[0] = ptr_new(f1)
    filtarray[1] = ptr_new(f2)

    frefarray = [1,0]
    
    nr = n_elements(filtarray)

    rpath = filepath('',root_dir=getenv('SIRTFZ_DIR'),subdirectory='results')
    openw, lun, rpath+'zmonte.log', /get_lun
;   printf, lun, [' ','Catalog','z = [0,1]','z = [1,2]','z > 2'], format='(A20,A21,A34,A33,A30)'
;   printf, lun, strjoin(replicate('-',155))

    ncatalog = 200
    mcnum = 100
    zmax = 6.0
    alpha = 2.0
    
    for i = 0L, nr-1L do begin

       filters = *filtarray[i]
       fref = frefarray[i]

       mc_catalog, filters, fref, zstats_all, zstats_zbins, ncatalog=ncatalog, $
         mcnum=mcnum, zmax=zmax, alpha=alpha, /postscript

       printf, lun, strjoin([filters,':',strcompress(string([zstats_all,zstats_zbins]),/remove)],' ')

    endfor
    free_lun, lun
    
stop
        
return
end
