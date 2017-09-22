pro dr2_densities

    ss = mrdfits('stripe82-dr12-stars.fits.gz',1)
    keep = where(ss.ra gt 320)
    ss = ss[keep]

    jj = mrdfits('decals-dr2-stripe82-dr12-stars.fits.gz',1,rows=keep)
    ww = where(strtrim(jj.brickname,2) ne '')
    help, ww, jj
    jj = jj[ww]
    ss = ss[ww]

    decals_to_maggies, jj, mm, ii, /decam_grz
    mm = mm[0:2,*] & ii = ii[0:2,*]
    good = where(mm[1,*] gt 0)
    help, good, jj

    sdss = ss[good]
    type = jj[good].type
    mag = 22.5-2.5*alog10(mm[*,good])

    psf = where(strtrim(type,2) eq 'PSF')
    


stop

return
end
