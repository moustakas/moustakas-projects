pro ndwfs_request_cutouts
; jm05jan04uofa
; jm10jan29ucsd - updated    
    
    agespath = ages_path(/analysis)
    outpath = ages_path(/thumbs)
    
    agesfile = agespath+'catalog.cat.noguidestars.fits.gz'
    splog, 'Reading '+agesfile
    weight = mrdfits(agespath+'catalog.spectweight.fits.gz',1)
    keep = where(weight.redshift gt 0.0)
    ages = mrdfits(agesfile,1,rows=keep)
    ngal = n_elements(ages)

    ages_id = string(keep,format='(I5.5)')
    
    ra = strtrim(im_dec2hms(ages.ra,/col),2)
    dec = repstr(strtrim(im_dec2hms(ages.dec,/col),2),'+','')
    rawidth = '0.5'

    band = 'I'
;   band = ['Bw','R','I','K']
    nband = n_elements(band)
    for iband = 0, nband-1 do begin
       if (file_test(outpath+band[iband],/dir) eq 0) then $
         spawn, 'mkdir -p '+outpath+band[iband], /sh
       fitsname = outpath+band[iband]+'/ages_'+ages_id+'.fits'

       scriptname = outpath+'cut_'+band[iband]+'.wget'
       openw, lun, scriptname, /get_lun
;      for ii = 4122L, ngal-1L do begin
       for ii = 0L, ngal-1L do begin
          printf, lun, 'wget "http://archive.noao.edu/ndwfs/cutout.php?ra='+$
            ra[ii]+'&dec='+dec[ii]+'&rawidth='+rawidth+'&decwidth=INDEF&fcsystem=J2000&filters='+band[iband]+'" '+$
            '-O "'+fitsname[ii]+'"'
       endfor
       free_lun, lun
       spawn, 'chmod +x '+scriptname, /sh
    endfor 
    
return
end
    
