pro sdss_emline_matrix, ll
; jm09feb23nyu - write out a matrix of emission-line properties
; for Hogg's archetype project

    inpath = sdss_path(/mpa_dr7)
    outpath = sdss_path(/proj)+'archetype/'

;   ll = read_sdss_vagc_mpa(sample='dr7',poststr='32',/ispec)
    if (n_elements(ll) eq 0L) then $
      ll = mrdfits(inpath+'gal_line_dr7_v5_2.fit.gz',1)
    
    matrix = [$
      [ll.sigma_balmer],$
      [ll.sigma_forbidden],$
      [ll.oii_flux],$
      [ll.oiii_5007_flux],$
      [ll.nii_6584_flux],$
      [ll.oi_6300_flux],$
      [ll.sii_6717_flux],$
      [ll.sii_6731_flux],$
      [ll.h_gamma_flux],$
      [ll.h_beta_flux],$
      [ll.h_alpha_flux]]

    errmatrix = [$
      [ll.sigma_balmer_err],$
      [ll.sigma_forbidden_err],$
      [ll.oii_flux_err],$
      [ll.oiii_5007_flux_err],$
      [ll.nii_6584_flux_err],$
      [ll.oi_6300_flux_err],$
      [ll.sii_6717_flux_err],$
      [ll.sii_6731_flux_err],$
      [ll.h_gamma_flux_err],$
      [ll.h_beta_flux_err],$
      [ll.h_alpha_flux_err]]

    linewave = [3727.42,5006.843,6584.04,6300.304,$
      6716.14,6730.81,4340.464,4861.325,6562.8]
    kl = k_lambda(linewave,/odonnel)

    mwrfits, matrix, outpath+'sdss_emline_matrix.fits', /create
    mwrfits, errmatrix, outpath+'sdss_emline_matrix.fits'
    mwrfits, [[linewave],[kl]], outpath+'sdss_emline_matrix.fits'
    spawn, 'gzip -f '+outpath+'sdss_emline_matrix.fits', /sh
    stop
return
end
