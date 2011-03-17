function init_sfrs_table, ngal
; initialize the data table
    
    sfrs = {$
      ages_id:        0L,$
      z:             0.0,$
      mips:            0,$ ; MIPS source?
      agn:             0,$ ; AGN?
      l1500:         0.0,$
      l1500_err:     0.0,$
      l1500_cor:     0.0,$ ; attenuation-corrected
      l1500_cor_err: 0.0,$ ; attenuation-corrected
      lir:           0.0,$ ; [erg/s]
      lir_err:       0.0,$ 
      l24:           0.0,$ ; [erg/s]
      l24_err:       0.0,$
      zmax_24:       0.0,$ ; redshift at which f24<f24_lim
      irx:           0.0,$ ; *linear* IRX from Meurer+99
      irx_err:       0.0,$
      a1500:         0.0,$ ; A(1500) from the flux-ratio method
      a1500_err:     0.0}
    sfrs = replicate(sfrs,ngal)

return, sfrs
end

pro build_cosmicimf_sfrs, sfrs
; jm10mar16ucsd - build the SFR estimates for the AGES/COSMICIMF
; sample

    cosmicimfpath = ages_path(/projects)+'cosmicimf/'

    splog, '#########################'
    splog, 'Building the AGES SFRs'
    sample = read_cosmicimf_sample()
    witt = read_cosmicimf_sample(/witt)
    ngal = n_elements(sample)

    sfrs = init_sfrs_table(ngal)
    sfrs.ages_id = sample.ages_id
    sfrs.z = sample.z

; identify AGN
    sfrs.agn = sample.x_match or ages_irac_agn(sample)

; UV luminosities [erg/s]
    sfrs.l1500 = alog10(1500D*sample.k_uvflux[0]*$
      4.0*!dpi*dluminosity(sample.z,/cm)^2/3.826D33)

; compute L(IR) and L(24) using the observed 24-micron flux densities
; and the average of the Chary & Elbaz (2001) and Dale & Helou (2002)
; infrared templates
    mips = where(sample.phot_mips24 gt 0.27,nmips)
    sfrs[mips].mips = 1
    f24 = sample[mips].phot_mips24
    f24_err = sample[mips].phot_mips24_err
    f24_limit = 0.27 ; limiting flux [mJy]

;    lir_r09 = im_f24tolir(sample.z,f24,f24_err=f24_err,/rieke,$
;      err_lir=lir_r09_err,l24=l24_r09,err_l24=l24_r09_err,$
;      f24_limit=f24_limit,zmax=zmax_r09)
;    sfrs.zmax_24 = zmax_r09
;    sfrs.lir = alog10(lir_r09)
;    sfrs.l24 = alog10(l24_r09)
;    sfrs.lir_err = lir_r09_err/lir_r09/alog(10.0)
;    sfrs.l24_err = l24_r09_err/l24_r09/alog(10.0)

    lir_ce01 = im_f24tolir(sample[mips].z,f24,f24_err=f24_err,/chary,$
      err_lir=lir_ce01_err,l24=l24_ce01,err_l24=l24_ce01_err,$
      f24_limit=f24_limit,zmax=zmax_ce01)
    lir_dh02 = im_f24tolir(sample[mips].z,f24,f24_err=f24_err,/dale,$
      err_lir=lir_dh02_err,l24=l24_dh02,err_l24=l24_dh02_err,$
      f24_limit=f24_limit,zmax=zmax_dh02)
    sfrs[mips].zmax_24 = total([[zmax_ce01],[zmax_dh02]],2)/2.0
    
    sfrs[mips].lir = total([[lir_ce01],[lir_dh02]],2)/2.0
    sfrs[mips].l24 = total([[l24_ce01],[l24_dh02]],2)/2.0
    for ii = 0L, nmips-1 do sfrs[mips[ii]].lir_err = djsig([lir_ce01[ii],lir_dh02[ii]])
    for ii = 0L, nmips-1 do sfrs[mips[ii]].l24_err = djsig([l24_ce01[ii],l24_dh02[ii]])
    sfrs[mips].lir_err = sfrs[mips].lir_err/sfrs[mips].lir/alog(10.0)
    sfrs[mips].l24_err = sfrs[mips].l24_err/sfrs[mips].l24/alog(10.0)
    sfrs[mips].lir = alog10(sfrs[mips].lir)
    sfrs[mips].l24 = alog10(sfrs[mips].l24)
    
; compute the infrared excess and A(1500) from the Witt & Gordon dust
; models; also correct L(1500) for attenuation
    sfrs[mips].irx = 10D^(sfrs[mips].lir-sfrs.l1500) ; keep this *linear*!
    sfrs[mips].irx_err = alog(10.0)*sfrs[mips].lir_err*sfrs[mips].irx

    linterp, witt.irx, witt.a1500, sfrs[mips].irx, a1500, missing=0.0
    sfrs[mips].a1500 = a1500
    sfrs[mips].l1500_cor = sfrs[mips].l1500 + 0.4*a1500

; get the errors; right now L(1500)_cor_err is effectively L(IR)_err,
; but once we get an error on L(1500) then this Monte Carlo technique
; will be the right thing to do
    nmonte = 100
    l1500_cor_monte = fltarr(nmonte,nmips)
    for ii = 0L, nmips-1L do begin
       l1500_monte = sfrs[mips[ii]].l1500 + randomn(seed,nmonte)*sfrs[mips[ii]].l1500_err
       irx_monte = sfrs[mips[ii]].irx + randomn(seed,nmonte)*sfrs[mips[ii]].irx_err
       linterp, witt.irx, witt.a1500, irx_monte, a1500_monte, missing=0.0
       sfrs[mips[ii]].a1500_err = djsig(a1500_monte)
       sfrs[mips[ii]].l1500_cor_err = djsig(l1500_monte + 0.4*a1500_monte)
    endfor
    
    im_mwrfits, sfrs, cosmicimfpath+'ages_cosmicimf_sfrs.fits', /clobber
       
return
end
