pro ages_select_qsos, allagesdust, allagesancillary
; jm06feb22uofa - select QSO's from the ALLAGESDUST spectral fits

    if (n_elements(allagesdust) eq 0L) then allagesdust = read_ages(ancillary=allagesancillary)

    linepath = ages_path(/analysis)
    linefile = 'ages_line_amplitudes.fits.gz'
    lineamps = mrdfits(linepath+linefile,1,/silent)

    rmscut = 1.5
    snrcut = 3.0
    
    mgiinarrow = where((lineamps.mgii_amplitude ne -1.0) and (lineamps.mgii_sigma lt 500.0) and $
                       (lineamps.mgii_amplitude gt (lineamps.mgii_medcont+rmscut*lineamps.mgii_rmscont)))
    struct_print, lineamps[mgiinarrow]
    mgiinarrow_specfit = read_ages_specfit(allagesdust[mgiinarrow].galaxy)
    ages_display_spectrum, allagesdust[mgiinarrow], allagesancillary[mgiinarrow], $
      specfit=mgiinarrow_specfit, labeltype=4L, /right, maxwave=5000.0, minwave=4800, setyrange=4L

    mgiibroad = where((lineamps.mgii_amplitude ne -1.0) and (lineamps.mgii_sigma ge 500.0) and $
                      (lineamps.mgii_amplitude gt (lineamps.mgii_medcont+rmscut*lineamps.mgii_rmscont)))
    struct_print, lineamps[mgiibroad]
    mgiibroad_specfit = read_ages_specfit(allagesdust[mgiibroad].galaxy)
    ages_display_spectrum, allagesdust[mgiibroad], allagesancillary[mgiibroad], $
      specfit=mgiibroad_specfit, labeltype=4L, /right, maxwave=4900.0

    hbetabroad = where((lineamps.hbeta_amplitude ne -1.0) and (allagesdust.h_beta_ew[0] gt 5.0) and $$ ;(lineamps.hbeta_sigma ge 200.0) and $
                       (lineamps.hbeta_amplitude gt (lineamps.hbeta_medcont+rmscut*lineamps.hbeta_rmscont)))
    struct_print, lineamps[hbetabroad]
    hbetabroad_specfit = read_ages_specfit(allagesdust[hbetabroad].galaxy)
    ages_display_spectrum, allagesdust[hbetabroad], allagesancillary[hbetabroad], $
      specfit=hbetabroad_specfit, labeltype=4L, /right, maxwave=6600.0

; --------------

    keep = where((lineamps.hbeta_amplitude ne -1.0) and (allagesdust.h_beta_ew[0] gt 5.0) and $
                 (lineamps.hbeta_amplitude gt (lineamps.hbeta_medcont+rmscut*lineamps.hbeta_rmscont)) and $
                 (lineamps.oii_amplitude gt (lineamps.oii_medcont+rmscut*lineamps.oii_rmscont)) and $
                 (lineamps.oiii_amplitude gt (lineamps.oiii_medcont+rmscut*lineamps.oiii_rmscont)) and $
                 (lineamps.mgii_sigma lt 500.0))
    keepindx = where((allagesdust[keep].h_beta[0]/allagesdust[keep].h_beta[1] gt snrcut) and $
                (allagesdust[keep].oii_3727[0]/allagesdust[keep].oii_3727[1] gt snrcut) and $
                (allagesdust[keep].oiii_5007[0]/allagesdust[keep].oiii_5007[1] gt snrcut) and $
                (strtrim(lineamps[keep].bpt_class,2) ne 'AGN'))
    srt = reverse(sort(lineamps[keep[keepindx]].z))
;   keepindx = where(lineamps[keep].z gt 0.5)
;   struct_print, lineamps[keep[keepindx]]
;   niceprint, allagesdust[keep[keepindx[srt]]].galaxy, allagesdust[keep[keepindx[srt]]].z_obj
    help, where(lineamps[keep[keepindx]].x_lum gt 0.0)
    keep_specfit = read_ages_specfit(allagesdust[keep[keepindx[srt]]].galaxy)
    ages_display_spectrum, allagesdust[keep[keepindx[srt]]], allagesancillary[keep[keepindx[srt]]], $
      specfit=keep_specfit, labeltype=4L, setyrange=3L, /right, maxwave=6600.0

; --------------
; test the things that are high-amplitude but low S/N
    
    test = where((lineamps.hbeta_amplitude ne -1.0) and (allagesdust.h_beta_ew[0] gt 5.0) and $
                 (lineamps.hbeta_amplitude gt (lineamps.hbeta_medcont+rmscut*lineamps.hbeta_rmscont)) and $
                 (lineamps.oii_amplitude gt (lineamps.oii_medcont+rmscut*lineamps.oii_rmscont)) and $
                 (lineamps.oiii_amplitude gt (lineamps.oiii_medcont+rmscut*lineamps.oiii_rmscont)) and $
                 (lineamps.mgii_sigma lt 500.0))
    testindx = where((allagesdust[test].h_beta[0]/allagesdust[test].h_beta[1] le snrcut) or $
                (allagesdust[test].oii_3727[0]/allagesdust[test].oii_3727[1] le snrcut) or $
                (allagesdust[test].oiii_5007[0]/allagesdust[test].oiii_5007[1] le snrcut))
    srt = reverse(sort(lineamps[test[testindx]].z))
;   testindx = where(lineamps[test].z gt 0.5)
;   struct_print, lineamps[test[testindx]]
;   niceprint, allagesdust[test[testindx[srt]]].galaxy, allagesdust[test[testindx[srt]]].z_obj
    help, where(lineamps[test[testindx]].x_lum gt 0.0)
    test_specfit = read_ages_specfit(allagesdust[test[testindx[srt]]].galaxy)
    ages_display_spectrum, allagesdust[test[testindx[srt]]], allagesancillary[test[testindx[srt]]], $
      specfit=test_specfit, labeltype=4L, setyrange=3L, maxwave=6600.0

; BPT plot    

    xtitle = 'log ([N II] \lambda6584 / H\alpha)_{obs}' & ytitle = 'log ([O III] \lambda5007 / H\beta)_{obs}'
    xrange = [-2.3,0.7] & yrange = [-1.3,1.4]

    hbewcut = where(allagesdust.h_beta_ew[0] gt 3.0)
    lineratio, allagesdust[hbewcut], 'NII_6584', 'H_ALPHA', 'OIII_5007', 'H_BETA', $
      x, xerr, y, yerr, snrcut=3.0

;   narrow = where((strtrim(allagesdust.bpt_pure_nii_class,2) eq 'AGN') and (allagesdust.h_beta_sigma[0] lt 200.0))
;   lineratio, allagesdust[narrow], 'NII_6584', 'H_ALPHA', 'OIII_5007', 'H_BETA', $
;     xnarrow, xerrnarrow, ynarrow, yerrnarrow, snrcut=3.0

;   broad = where((strtrim(allagesdust.bpt_pure_nii_class,2) eq 'AGN') and (allagesdust.h_beta_sigma[0] ge 200.0))
    broad = where((allagesdust.h_beta_sigma[0] ge 200.0) and (allagesdust.h_beta_ew[0] gt 3.0))
    lineratio, allagesdust[broad], 'NII_6584', 'H_ALPHA', 'OIII_5007', 'H_BETA', $
      xbroad, xerrbroad, ybroad, yerrbroad, snrcut=3.0

    ages_lineplot, x, y, xerr, yerr, plottype=3, postscript=postscript, $
      xtitle=xtitle, ytitle=ytitle, xrange=xrange, yrange=yrange
;   ages_lineplot, xnarrow, ynarrow, xerrnarrow, yerrnarrow, /overplot, agescolor='blue'
    ages_lineplot, xbroad, ybroad, xerrbroad, yerrbroad, /overplot, agescolor='green'
    
;   agn = where((strtrim(allagesdust.bpt_pure_nii_class,2) eq 'AGN'))
;   plothist, allagesdust[agn].h_beta_sigma[0]
    agn = where((strtrim(allagesdust.bpt_pure_nii_class,2) eq 'AGN') and (lineamps.hbeta_amplitude ne -1.0) and $
                (lineamps.hbeta_sigma ge 200.0) and (lineamps.hbeta_amplitude gt (lineamps.hbeta_medcont+rmscut*lineamps.hbeta_rmscont)))
    struct_print, lineamps[agn]
    hbetabroad_specfit = read_ages_specfit(allagesdust[hbetabroad].galaxy)
    ages_display_spectrum, allagesdust[hbetabroad], allagesancillary[hbetabroad], $
      specfit=hbetabroad_specfit, labeltype=4L, /right, maxwave=6600.0


stop
    
    snrcut = 5.0
    
    loz = where(((allagesdust.h_beta[0]/allagesdust.h_beta[1] gt snrcut) and (allagesdust.z_obj lt 0.32) and $
      (allagesdust.h_beta_sigma[0] gt 150.0)) or ((strtrim(allagesancillary.spzbest_class,2) eq 'QSO') and $
      (allagesdust.z_obj lt 0.32)),nloz)
    niceprint, allagesdust[loz].galaxy, allagesdust[loz].z_obj, allagesdust[loz].h_beta_sigma[0], $
      allagesdust[loz].h_beta_ew[0], allagesancillary[loz].spzbest_class
    help, loz

; ---------------------------------------------------------------------------    
    
    midz = where(((allagesdust.mgii_2800[0]/allagesdust.mgii_2800[1] gt snrcut) and $
                  (allagesdust.h_beta[0]/allagesdust.h_beta[1] gt snrcut) and (allagesdust.z_obj ge 0.32) and $
                  (allagesdust.z_obj lt 0.75) and (allagesdust.h_beta_sigma[0] gt 150.0) and $
                  (allagesdust.mgii_2800_sigma[0] gt 500.0) and (allagesdust.h_beta_ew[0] gt 2.0) and $
                  (allagesdust.mgii_2800_ew[0] gt 0.0) and (allages),nmidz); or $
;                ((strtrim(allagesancillary.spzbest_class,2) eq 'QSO') and (allagesdust.z_obj ge 0.32) and $
;                 (allagesdust.z_obj lt 0.75)),nmidz)
    print, minmax(allagesdust[midz].z_obj)
    niceprint, allagesdust[midz].galaxy, long(allagesancillary[midz].fluxing_problem), allagesdust[midz].z_obj, $
      allagesdust[midz].h_beta_sigma[0], allagesdust[midz].oii_3727[0]/allagesdust[midz].oii_3727[1], $
      allagesdust[midz].mgii_2800_sigma[0], allagesdust[midz].h_beta_ew[0], $
      allagesdust[midz].mgii_2800_ew[0], allagesancillary[midz].spzbest_class
    help, midz

    indx = lindgen(nmidz)
    indx = where(allagesdust[midz].h_beta_ew[0] lt 5.0)
    indx = where(allagesdust[midz].oii_3727[0]/allagesdust[midz].oii_3727[1] lt snrcut)
    specfit = read_ages_specfit(allagesdust[midz[indx]].galaxy)
    ages_display_spectrum, allagesdust[midz[indx]], allagesancillary[midz[indx]], $
      specfit=specfit, labeltype=4L, /right, maxwave=5000.0
    
; ---------------------------------------------------------------------------    

    
    loqso = where(((allagesdust.mgii_2800[0]/allagesdust.mgii_2800[1] gt snrcut) and $
      (allagesdust.h_beta[0]/allagesdust.h_beta[1] gt snrcut) and (allagesdust.z_obj lt 0.75) and $
      ((allagesdust.mgii_2800_sigma[0] gt 500.0) or (allagesdust.h_beta_sigma[0] gt 300.0))))
      (strtrim(allagesancillary.spzbest_class,2) eq 'QSO'))
    niceprint, allagesdust[loqso].galaxy, allagesdust[loqso].z_obj, allagesdust[loqso].mgii_2800_sigma[0], $
      allagesdust[loqso].h_beta_sigma[0], allagesancillary[loqso].spzbest_class
    help, loqso

    lo_specfit = read_allagesdust_specfit(allagesdust[loqso[0:100]].galaxy)
    allagesdust_display_spectrum, allagesdust[loqso[0:100]], allagesancillary[loqso[0:100]], specfit=lo_specfit

    stop    
    
return
end
    
