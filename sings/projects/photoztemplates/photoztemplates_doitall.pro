pro photoztemplates_doitall
    unpack_photoztemplates_phot
    photoztemplates_kcorrect, /clobber
    photoztemplates_qaplot
    photoztemplates_qaplot, /emlines
    photoztemplates_write_kcorrect_models
    spawn, ['tar czvf photoztemplates.tar.gz kcorrect_models.fits.gz '+$
      'kcorrect_models_emlines.fits.gz kcorr.fits.gz kcorr_emlines.fits.gz '+$
      'qaplot_photoztemplates.ps.gz qaplot_photoztemplates_emlines.ps.gz'], /sh
return
end
