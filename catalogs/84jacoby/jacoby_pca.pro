pro jacoby_pca
; jm03jul31uofa

    path = filepath('',root_dir=getenv('CATALOGS_DIR'),subdirectory='jacoby_atlas/eigen')
    pushd, path
    flist = findfile('*.fits',count=fcount)
    for k = 0L, fcount-1L do begin

       objflux1 = readfits(flist[k],h,/silent)
       if k eq 0L then objflux = objflux1 else objflux = [ [objflux], [objflux1] ]
       
    endfor

    objlam = make_wave(h)
    objivar = objflux*0.0+1.0
    objloglam = alog10(objlam)

    zfit = 0.0
    
    pcaflux = pca_solve(objflux,objivar,objloglam,zfit,$
      wavemin=wavemin,wavemax=wavemax,newloglam=newloglam,$
      maxiter=maxiter,niter=niter,nkeep=nkeep,nreturn=nreturn,$
      eigenval=eigenval,acoeff=acoeff,outmask=outmask,$
      usemask=usemask,_extra=extra)

    wave = 10.0^newloglam
    crval1 = wave[0]
    cd1_1 = wave[1]-wave[0]
    eigenres = 4.5
    
    mkhdr, header, pcaflux
    sxaddpar, header, 'OBJECT', 'Jacoby PCA templates '+im_today()
    sxaddpar, h, 'EXTEND', 'N'
    sxaddpar, header, 'CRVAL1', crval1, ' central wavelength of first pixel'
    sxaddpar, header, 'CD1_1', cd1_1, ' dispersion [Angstrom/pixel]'
    sxaddpar, header, 'CRPIX1', 0, ' starting pixel (0-indexed)'
    sxaddpar, header, 'CTYPE1', 'LINEAR'
    sxaddpar, header, 'DC-FLAG', 0, ' log-linear flag'
    sxaddpar, header, 'EIGENRES', eigenres, 'resolution [Angstrom]'

    if keyword_set(write) then begin

       tpath = filepath('',root_dir=getenv('ISPEC_DIR'),subdirectory='templates')
       tfile = 'jacoby_pca.fits'
       print, 'Writing '+tpath+tfile+'.'
       mwrfits, pcaflux, tpath+tfile, header, /create

    endif

return
end    
