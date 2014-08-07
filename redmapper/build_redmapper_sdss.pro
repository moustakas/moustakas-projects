pro build_redmapper_sdss
; jm13mar28siena - build the SDSS catalog for the sample using
; CasJobs; be sure to write the output file from CasJobs as
; 'redmapper_'+ver+'_sdss.fits' 

    path = redmapper_path(version=ver)
    cat = mrdfits(path+'dr8_run_redmapper_'+ver+$
      '_lgt5_catalog_members.fits.gz',1)
    ngal = n_elements(cat)

    out = struct_addtags(replicate({redmapper_id: 0L, casid: long64(0)},ngal),cat)
    out.redmapper_id = lindgen(ngal)
    out.casid = photoid2casid(cat.photoid)

    nchunk = 4
    chunksize = long(ngal/float(nchunk))

    for ichunk = 0, nchunk-1 do begin
       i1 = ichunk*chunksize
       i2 = (((ichunk*chunksize+chunksize)<ngal)-1L)>0L
       nindx = i2-i1+1L 
       indx = lindgen(nindx)+i1
     
       outfile = '~/tmp/redmapper_'+ver+'_casjobs_'+string(ichunk,format='(I2.2)')+'.dat'
       splog, 'Writing '+outfile
       openw, lun, outfile, /get_lun
;      printf, lun, '# redmapper_id casid ra dec'
;      struct_print, struct_trimtags(out,select=['redmapper_id','casid','ra','dec']), $
       printf, lun, '# redmapper_id casid'
       struct_print, struct_trimtags(out[indx],select=['redmapper_id','casid']), $
         lun=lun, ddigit=12, /no_head
       free_lun, lun
    endfor

return
end
    
