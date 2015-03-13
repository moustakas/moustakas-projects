pro merge_tractor_dr1, merge=merge, specz=specz
; jm15mar11siena - merge the DR1 tractor catalogs    

    dr1_dir = getenv('DECALS_DIR')+'/'
    isedfit_dir = getenv('DECALS_DIR')+'/isedfit/dr1/'

    if keyword_set(merge) then begin
       allbrick = file_basename(file_search(dr1_dir+'tractor/*',/test_dir,count=nbrick))
       for ii = 0L, nbrick-1 do begin
          catfile = file_search(dr1_dir+'tractor/'+allbrick[ii]+'/tractor-*.fits',count=ncat)
          for ic = 0L, ncat-1 do begin
             cat1 = mrdfits(catfile[ic],1)
             if n_elements(cat) eq 0L then cat = cat1 else cat = [cat,cat1]
          endfor
       endfor
       keep = where(cat.brick_primary eq 'T')
       im_mwrfits, cat[keep], dr1_dir+'tractor_dr1.fits', /clobber
    endif

    if keyword_set(specz) then begin
       cat = mrdfits(dr1_dir+'tractor_dr1.fits.gz',1)
       specz = mrdfits(getenv('
       
    endif
    
stop

return
end
