pro compare_sdss_deep2groth
; jm14feb28siena - compare [OII] fluxes of objects in the DEEP2 Groth
; Strip (Field 1) to those measured from some special SDSS plates

    outpath = getenv('IM_PROJECTS_DIR')+'/desi/targeting/'

    deep2 = read_deep2_zcat()
    ppxf = read_deep2(/ppxf)
    kcorr = read_deep2(/kcorr)

    sdss = mrdfits(outpath+'sdss-plug-grothstrip.fits.gz',1)
    line = mrdfits(outpath+'sdss-zline-grothstrip.fits.gz',1)
    line = reform(line,31,n_elements(line)/31)
    
    zans = mrdfits(outpath+'sdss-zans-grothstrip.fits.gz',1)
    spherematch, deep2.ra, deep2.dec, sdss.ra, sdss.dec, 1D/3600.0, m1, m2
    keep = where(zans[m2].zwarning eq 0 and (line[6,*].linearea_err gt 0 or $
      line[7,*].linearea_err gt 0) and (ppxf[m1].oii_3727_1[1] gt 0 or $
      ppxf[m1].oii_3727_2[1] gt 0),ngal)
    help, ngal
    m1 = m1[keep]
    m2 = m2[keep]

    deep2 = deep2[m1]
    ppxf = ppxf[m1]
    kcorr = kcorr[m1]

    sdss = sdss[m2]
    zans = zans[m2]
    line = line[*,m2]

    soii = reform(line[6,*].linearea+line[7,*].linearea)
    soiierr = reform(sqrt(line[6,*].linearea_err^2+line[7,*].linearea_err^2))

    doii = ppxf.oii_3727_1[0]+ppxf.oii_3727_2[0]
    doiierr = sqrt(ppxf.oii_3727_1[1]^2+ppxf.oii_3727_2[1]^2)

stop    

return
end



