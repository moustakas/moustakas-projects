pro build_wisesfrs_sample, clobber=clobber
; jm11apr19ucsd - build the final science sample for the WISE SFR
; calibration paper 
    
    sfrspath = wise_path(/wisesfrs)
    info = mrdfits(sfrspath+'wisesfrs_info.fits.gz',1)
    wise = mrdfits(sfrspath+'wisesfrs_wisecat.fits.gz',1)

; get L(IR)    
    wise_to_maggies, wise, maggies, ivarmaggies, /nodust
    lir_r09 = im_wise2lir(info.zdist,maggies[2:3,*],ivarmaggies[2:3,*],err_lir=lirerr_ce01,/rieke,/zlog,/debug)

stop    

return
end
    
