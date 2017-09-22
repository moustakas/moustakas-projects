pro choose_dr3a

    jj = mrdfits('~/repos/git/legacysurvey/legacypipe-dir/decals-bricks.fits',1)

; DEEP2/Field 4
    ff = mrdfits('~/research/projects/desi/target/analysis/truth/deep2-field4.fits.gz',1)
    ra = minmax(ff.ra_deep)+0.25/2*[-1,1]
    dec = minmax(ff.dec_deep)+0.25/2*[-1,1]

    djs_plot, ff.ra_deep, ff.dec_deep, xsty=1, ysty=1, $
      xrange=ra, yrange=dec, $
      psym=6, symsize=2
    djs_oplot, jj.ra, jj.dec, psym=6, color='orange'

    ww = where(jj.ra gt ra[0] and jj.ra lt ra[1] and $
      jj.dec gt dec[0] and jj.dec lt dec[1])
    niceprint, jj[ww].brickname
    
; VIPERS/W1
    ff = mrdfits('~/research/projects/desi/target/analysis/truth/vipers-w1.fits.gz',1)
    ra = minmax(ff.alpha)+0.25/2*[-1,1]
    dec = minmax(ff.delta)+0.25/2*[-1,1]

    djs_plot, ff.alpha, ff.delta, xsty=1, ysty=1, $
      xrange=ra, yrange=dec, $
      psym=6, symsize=2
    djs_oplot, jj.ra, jj.dec, psym=6, color='orange'

    ww = where(jj.ra gt ra[0] and jj.ra lt ra[1] and $
      jj.dec gt dec[0] and jj.dec lt dec[1])
    niceprint, jj[ww].brickname
    


return
end

