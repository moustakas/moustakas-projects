pro build_redmapper_unwise, query=query, parse=parse
; jm13mar28siena - match the REDMAPPER/v5.10 catalog to
; Dustin's unWISE catalogs

    ver = 'v5.10'
    path = '~/'
;   path = '/global/u2/i/ioannis/'
;   path = redmapper_path(version=ver)
;   unwisepath = '/global/project/projectdirs/desi/users/dstn/sdssv4-pobj/301/'
    unwisepath = '/moustakas-archive/sdssv4-pobj/301/'

    ngal = 11235932L
    chunksize = 1000000L
    nchunk = ceil(ngal/float(chunksize))

;   cat = mrdfits(path+'redmapper_'+ver+'_photoid.fits.gz',1)
;   ngal = n_elements(cat)
;   allfile = strtrim(cat.file,2)

; build the output structure
    out_template = mrdfits(unwisepath+'1000/1/photoWiseForced-001000-1-0027.fits',row=0,1)
    
    for ichunk = 10, nchunk-1 do begin
;   for ichunk = 0, nchunk-1 do begin
       i1 = ichunk*chunksize
       i2 = (((ichunk*chunksize+chunksize)<ngal)-1L)>0L
       nindx = i2-i1+1L 
       indx = lindgen(nindx)+i1
       
       cat = mrdfits(path+'redmapper_'+ver+'_photoid.fits.gz',1,rows=indx)
       allfile = strtrim(cat.file,2)
       ncat = n_elements(cat)
       
       out = replicate(out_template,ncat)
       
       ufile = allfile[uniq(allfile,sort(allfile))]
       nfile = n_elements(ufile)
       splog, 'Number of unique files: ', nfile
       
       for ii = 0L, nfile-1 do begin
          if (ii mod 100) eq 0 then print, ichunk, ii
          these = where(ufile[ii] eq allfile,nthese)
          unwise1 = mrdfits(unwisepath+ufile[ii],1,row=cat[these].id-1,/silent)
          
          out[these] = struct_trimtags(unwise1,except=['w3_*','w4_*'])
       endfor
       
       mwrfits, out, path+'redmapper_'+ver+'_unwise_'+string(ichunk,format='(I2.2)')+'.fits', /create
    endfor

stop    
    
    ss = mrdfits('redmapper_v5.10_unwise_00.fits',1)
    ff = file_search('redmapper_*.fits',count=cc)
    for ii=1,cc-1 do ss = [temporary(ss),mrdfits(ff[ii],1)]
    im_mwrfits, ss, path+'redmapper_v5.10_unwise.fits', /clobber

return
end    
