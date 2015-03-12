pro merge_tractor_dr1

    dr1path = '/Users/ioannis/research/projects/decals/dr1/'
    allbrick = file_basename(file_search(dr1path+'tractor/*',/test_dir,count=nbrick))
    for ii = 0L, nbrick-1 do begin
       catfile = file_search(dr1path+'tractor/'+allbrick[ii]+'/tractor-*.fits',count=ncat)
       for ic = 0L, ncat-1 do begin
          cat1 = mrdfits(catfile[ic],1)
          if n_elements(cat) eq 0L then cat = cat1 else cat = [cat,cat1]
       endfor
    endfor

    ww = where(cat.brick_primary eq 'T')
    mwrfits, cat[ww], 'tractor_dr1.fits', /create

stop

return
end
