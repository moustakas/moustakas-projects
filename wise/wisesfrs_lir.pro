pro wisesfrs_lir, clobber=clobber, debug=debug
; jm11apr19ucsd - compute L(IR) for the parent sample
    
    sfrspath = wise_path(/wisesfrs)
    atlas = mrdfits(sfrspath+'wisesfrs_atlas.fits.gz',1)
    wise = mrdfits(sfrspath+'wisesfrs_wisecat.fits.gz',1)
    kcorr = mrdfits(sfrspath+'wisesfrs_kcorr.fits.gz',1)
    ngal = n_elements(atlas)

    out = struct_trimtags(atlas,select=['ra','dec','zdist'])
    out = struct_addtags(out,replicate({$
      lir_ce01: 0.0, lir_ce01_err: 0.0, chi2_ce01: 0.0, modelindx_ce01: 0.0,$
      lir_dh02: 0.0, lir_dh02_err: 0.0, chi2_dh02: 0.0, modelindx_dh02: 0.0,$
      fmnuv: 0.0, beta: 0.0, l1500: 0.0, irx_ce01: 0.0, irx_dh02: 0.0},ngal))
;     lir_r09: 0.0, lir_r09_err: 0.0},ngal))
    
    wise_to_maggies, wise, maggies, ivarmaggies

; Chary & Elbaz 01
    lir = im_wise2lir(atlas.zdist,maggies,ivarmaggies,$
      err_lir=lirerr,/chary,/zlog,chi2=chi2,debug=debug,$
      model_indx=indx)
    out.lir_ce01 = lir
    out.chi2_ce01 = chi2
    out.modelindx_ce01 = indx

; Dale & Helou 02
    lir = im_wise2lir(atlas.zdist,maggies,ivarmaggies,$
      err_lir=lirerr,/dale,/zlog,chi2=chi2,debug=debug,$
      model_indx=indx)
    out.lir_dh02 = lir
    out.chi2_dh02 = chi2
    out.modelindx_dh02 = indx

;; Rieke+09
;    lir = im_wise2lir(atlas.zdist,maggies,ivarmaggies,$
;      err_lir=lirerr,/rieke,/zlog,chi2=chi2,debug=debug)
;    out.lir_r09 = lir
;    out.chi2_r09 = chi2

; compute L(1500), beta, IRX, etc.    
    out.fmnuv = kcorr.k_galex_absmag_00[0]-kcorr.k_galex_absmag_00[1]
    out.beta = 2.32*out.fmnuv-2

    out.l1500 = 1500D*kcorr.k_uvflux[0]*$
      4.0*!dpi*dluminosity(atlas.zdist,/cm)^2/3.826D33
    out.irx_ce01 = alog10(out.lir_ce01/out.l1500)
    out.irx_dh02 = alog10(out.lir_dh02/out.l1500)
    
    outfile = sfrspath+'wisesfrs_lir.fits'
    im_mwrfits, out, outfile, clobber=clobber
    
stop    

return
end
    
