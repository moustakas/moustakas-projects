pro build_icl_mosaics
; jm13may23siena - build the mosaics for each cluster

    redmapper_path = getenv('HOME')+'/redmapper/'
    outpath = redmapper_path+'decamlegacy/'
;   outpath = redmapper_path+'mosaics/'
    
; global parameters
    rerun = 301
    npixround = 8
    pixscale = 0.4D/3600D
    rad = 2D ; [Mpc]
    radeg = 206265D/dangular(0.33,/Mpc)/3600D*rad ; [degrees]

; read the cluster catalog and choose the best    
;   allcat = mrdfits(redmapper_path+'dr8_run_redmapper_v5.2_lgt20_catalog.fits.gz',1)
;   these = where(allcat.z lt 0.33 and allcat.lambda_chisq gt 120,nthese)
;   cat = allcat[these]
    cat = mrdfits(redmapper_path+'redmapper_decam_14aug.fits.gz',1)

    pushd, outpath
    t0 = systime(1)
    for ii = 0L, 0 do begin
;   for ii = 0L, nthese-1 do begin
       prefix = 'cl_'+string(cat[ii].mem_match_id,format='(I5.5)')
       jpgfile = outpath+prefix+'.jpg'
       smosaic_make_jpg, cat[ii].ra, cat[ii].dec, radeg, radeg, jpgfile, $
         rerun=rerun, objlist={run:0, camcol:0, field:0, id:0, rerun:''}, $
         /globalsky, /ivarout, prefix=prefix, pixscale=pixscale
stop
;      smosaic_make, cat[ii].ra, cat[ii].dec, radeg, radeg, rerun=rerun, $
;        objlist={run:0, camcol:0, field:0, id:0, rerun:''}, $
;        /globalsky, /ivarout, prefix=prefix, pixscale=pixscale
    endfor
    splog, 'Total time = ', (systime(1)-t0)/60.0
    popd

return
end
    
