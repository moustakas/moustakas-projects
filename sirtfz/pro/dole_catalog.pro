pro dole_catalog, ncatalog=ncatalog
; jm01sep10uofa

    if not keyword_set(ncatalog) then ncatalog = 100L
    
    dolepath = '/home/ioannis/sirtf/simulations/savesets/'

    print, 'Restoring '+dolepath+'lfz.save.'
    restore, dolepath+'lfz.save'
;   print, 'Restoring '+dolepath+'multiwavelength_sources_1024.save.'
;   restore, dolepath+'multiwavelength_sources_1024.save'

    filters = ['IRAC 3.6','IRAC 4.5','IRAC 5.8','IRAC 8','MIPS 24','MIPS 70','MIPS 160']
    findx = [0,1,2,3,4,5,6] ; sirtf filters
    nbands = n_elements(findx)

    lfzsize = size(lfz,/dimension)
    nseds = lfzsize[0]
    nz = lfzsize[2]

    flimit = [5.0E-4,7.2E-4,2.3E-3,3.1E-3,0.03,0.2,6.0] ; limiting flux [mJy]
    lambda0 = wavelengths_lfz[findx]*1E6 ; [micron]
    
    path = filepath('',root_dir=getenv('SIRTFZ_DIR'),subdirectory='catalogs')
    file = 'dole.cat'
    print, 'Opening '+path+file+'...'
    
    openw, lun, path+file, /get_lun
    printf, lun, '# H. Dole simulated catalog: ncat = '+strn(ncatalog)
    header = ['# '+filters[0]+',',filters[1:nbands-2L]+',',filters[nbands-1L]]
    printf, lun, header

    nsources = 0L
    while nsources le ncatalog do begin

       zindx = floor(randomu(seed,1)*nz)    ; random redshift
       sindx = floor(randomu(seed,1)*nseds) ; random SED type (L_IR)

;      z = z_lfz[zindx]
       z = z_lfz[673]

       flux = (reform(lfz[sindx,*,zindx]))[findx]*1E3 ; [mJy]
       ferror = sqrt(flimit^2 + (0.05*flux)^2)

       dflux = flux + randomn(seed,nbands)*ferror ; Gaussian-perturb the flux

       ploterror, lambda0, dflux, ferror, ps=4, xsty=3, ysty=3
       
       print, format='("Writing source ",I0," of ",I0,".",A1,$)', nsources, ncatalog, string(13b)
       printf, lun, dflux, ferror, z, format='('+strn(nbands)+'E15.5,'+strn(nbands)+'E15.5,3x,D0.0)'

       stop
       
       nsources = nsources + 1L

    endwhile
    print
    free_lun, lun
    
return
end
