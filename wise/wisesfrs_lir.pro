pro wisesfrs_lir, clobber=clobber, debug=debug
; jm11apr19ucsd - compute L(IR) for the parent sample
    
    sfrspath = wise_path(/wisesfrs)
    atlas = mrdfits(sfrspath+'wisesfrs_atlas.fits.gz',1)
    wise = mrdfits(sfrspath+'wisesfrs_wisecat.fits.gz',1)
    ngal = n_elements(atlas)

    out = struct_trimtags(atlas,select=['ra','dec','zdist'])
    out = struct_addtags(out,replicate({lir_ce01: 0.0, lir_ce01_err: 0.0, $
      lir_dh02: 0.0, lir_dh02_err: 0.0, lir_r09: 0.0, lir_r09_err: 0.0},ngal))
    
    wise_to_maggies, wise, maggies, ivarmaggies

; Chary & Elbaz 01
    lir = im_wise2lir(atlas.zdist,maggies[2:3,*],ivarmaggies[2:3,*],$
      err_lir=lirerr,/chary,/zlog,debug=debug)
    out.lir_ce01 = lir

; Dale & Helou 02
    lir = im_wise2lir(atlas.zdist,maggies[2:3,*],ivarmaggies[2:3,*],$
      err_lir=lirerr,/dale,/zlog,debug=debug)
    out.lir_dh02 = lir

; Rieke+09
    lir = im_wise2lir(atlas.zdist,maggies[2:3,*],ivarmaggies[2:3,*],$
      err_lir=lirerr,/rieke,/zlog,debug=debug)
    out.lir_r09 = lir

    outfile = sfrspath+'wisesfrs_lir.fits'
    im_mwrfits, out, outfile, clobber=clobber
    
stop    

return
end
    
