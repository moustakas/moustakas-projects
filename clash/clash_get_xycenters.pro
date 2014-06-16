pro clash_get_xycenters
; jm14jun14siena - get the xy coordinates corresponding to the cluster
;   centers in CLASH_SAMPLE.SEX for Adi

    clash = rsex(getenv('CLASH_DIR')+'/clash_sample.sex')
    ncl = n_elements(clash)

    out = struct_addtags(im_struct_trimtags(clash,select=['shortname','ra','dec'],$
      newtags=['cluster','ra','dec']),replicate({xcen: 0D, ycen: 0D},ncl))
    
    for ic = 0, ncl-1 do begin
       cluster = strtrim(clash[ic].shortname,2)
       splog, cluster
       mosaicpath = getenv('CLASH_ARCHIVE')+'/'+strtrim(clash[ic].dirname,2)+$
         '/HST/images/mosaicdrizzle_image_pipeline/scale_65mas/'
       imfile = file_search(mosaicpath+cluster+'_mosaic_065mas_wfc3ir_f160w_drz_????????.fits*')
       if file_test(imfile) eq 0 then stop
       hdr = headfits(imfile)
       adxy, hdr, 15D*hms2dec(clash[ic].ra), hms2dec(clash[ic].dec), xx, yy
       out[ic].xcen = xx
       out[ic].ycen = yy
    endfor

    struct_print, out
    struct_print, out, file='clash_xycenters.txt'

return
end

