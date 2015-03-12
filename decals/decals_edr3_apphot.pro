pro decals_edr3_apphot
; jm15feb25siena - do aperture photometry on the residual images
; (centered on each source) in each brick

    edrpath = '/Users/ioannis/research/projects/decals/edr3/'
;   edrpath = '/global/work/decam/cats/edr3/'
    outpath = edrpath+'apphot/'

    nrad = 15
    radius = findgen(nrad)+1 ; pixels

    band = ['g','r','z']
    nband = n_elements(band)

    out1 = {$
      apphot:            fltarr(3,nrad),$
      apphot_ivar:       fltarr(3,nrad),$
      apphot_resid:      fltarr(3,nrad),$
      apphot_resid_ivar: fltarr(3,nrad)}
    
    allbrick = file_basename(file_search(edrpath+'tractor/*',/test_dir,count=nbrick))
    for ii = 0L, nbrick-1 do begin
       catfile = file_search(edrpath+'tractor/'+allbrick[ii]+'/tractor-*.fits',count=ncat)
       outcatfile = repstr(catfile,'tractor-','tractor2-')
       for ic = 0L, ncat-1 do begin
          cat = mrdfits(catfile[ic],1)
          nobj = n_elements(cat)
          outcat = struct_addtags(cat,replicate(out1,nobj))

          brick = strtrim(cat[0].brickname,2)
          coaddpath = edrpath+'coadd/'+allbrick[ii]+'/'+brick+'/'

          for ib = 0, nband-1 do begin
             imfile = coaddpath+'decals-'+brick+'-image-'+band[ib]+'.fits'
             image = mrdfits(imfile,0,hdr)
             invvar = mrdfits(repstr(imfile,'-image','-invvar'),0)
             model = mrdfits(repstr(imfile,'-image','-model'),0)
             resid = image - model

             flux = djs_phot(cat.bx,cat.by,radius,skyrad,image,$
               invvar,calg='none',flerr=ferr,peakval=peak)

stop

             flux_resid = djs_phot(cat.bx,cat.by,radius,skyrad,resid,$
               invvar,calg='none',flerr=ferr_resid,peakval=peak_resid)

             outcat.apphot[ib,*] = reform(transpose(flux),1,nrad,nobj)
             outcat.apphot_resid[ib,*] = reform(transpose(flux_resid),1,nrad,nobj)

             outcat.apphot_ivar[ib,*] = reform(transpose(1.0/ferr^2*(ferr gt 0)),1,nrad,nobj)
             outcat.apphot_resid_ivar[ib,*] = reform(transpose(1.0/ferr_resid^2*(ferr_resid gt 0)),1,nrad,nobj)
          endfor 
          mwrfits, outcat, outcatfile[ic], /create
       endfor 
    endfor

return
end
