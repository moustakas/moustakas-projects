pro photoztemplates_doitall
    unpack_photoztemplates_phot
    photoztemplates_kcorrect, /clobber
    photoztemplates_qaplot
    photoztemplates_write_kcorrect_models
    spawn, ['tar czvf photoztemplates.tar.gz kcorrect_models.fits.gz '+$
      'kcorr.fits.gz qaplot_photoztemplates.ps.gz'], /sh
return
end
