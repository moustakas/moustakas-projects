pro merge_tractor_dr1, merge=merge, specz=specz
; jm15mar11siena - merge the DR1 tractor catalogs    

    dr1_dir = getenv('DECALS_DIR')+'/'
    isedfit_dir = getenv('IM_ARCHIVE_DIR')+'/decam/isedfit/dr1/'

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
       im_mwrfits, cat[keep], isedfit_dir+'tractor-dr1.fits', /clobber
    endif

; match against the specz catalog    
    if keyword_set(specz) then begin
       cat = mrdfits(isedfit_dir+'tractor_dr1.fits.gz',1)
       specz = mrdfits(getenv('IM_ARCHIVE_DIR')+'/data/sdss/dr12/specObj-dr12.fits',1)
       spherematch, specz.plug_ra, specz.plug_dec, cat.ra, cat.dec, 1D/3600, m1, m2

       keep = where(strtrim(specz[m1].class,2) eq 'GALAXY' and specz[m1].zwarning eq 0 and $
         specz[m1].z ge 0.05 and specz[m1].z le 0.7,ngal)
       out = struct_addtags(cat[m2[keep]],replicate({z: 0.0},ngal))
       out.z = specz[m1[keep]].z
       
       im_mwrfits, out, isedfit_dir+'tractor-dr1-specz.fits', /clobber
       
stop       
       
    endif
    
stop

return
end
