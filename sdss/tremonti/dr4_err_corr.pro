
function dr4_err_corr, sl

; From Jarle (for dr4_v51b)

sl.OII_3726_flux_err = sl.OII_3726_flux_err * 1.7346
sl.OII_3729_flux_err = sl.OII_3729_flux_err * 1.7346
sl.NeIII_3869_flux_err = sl.NeIII_3869_flux_err * 1.2456
sl.H_beta_flux_err = sl.H_beta_flux_err * 1.8458
sl.OIII_4363_flux_err = sl.OIII_4363_flux_err * 1.1730
sl.OIII_4959_flux_err = sl.OIII_4959_flux_err * 1.3287
sl.OIII_5007_flux_err = sl.OIII_5007_flux_err * 1.4435
sl.HeI_5876_flux_err = sl.HeI_5876_flux_err * 1.2744
sl.OI_6300_flux_err = sl.OI_6300_flux_err * 1.1837
sl.H_alpha_flux_err = sl.H_alpha_flux_err * 2.5679
sl.NII_6584_flux_err = sl.NII_6584_flux_err * 1.8885
sl.SII_6717_flux_err = sl.SII_6717_flux_err * 1.6323
sl.SII_6731_flux_err = sl.SII_6731_flux_err * 1.6323

return, sl

end
