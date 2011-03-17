pro k04_test, k = k, s = s, ns=ns

    dir = '/home/tremonti/spectra/john/sfr/'
    dir = './'

    if not keyword_set(k) or not keyword_set(s) then begin
       k = rsex(dir + '04kewley_table1.txt')
; s = mrdfits(dir + 'specfit/nfgs_int_speclinefit.fits.gz', 1)
       s = read_nfgs(nfgsnodust=ns)
       s = mrdfits(dir + 'specfit/nfgs_int_speclinefit.fits.gz', 1)
    endif

    nk = n_elements(k)

;----------------------------------
; Matching

    nfgs_id = long(strmid(s.drift_file, 0, 3))
    sok = cmset_op(nfgs_id, 'and', k.id, /index)

;-----------------------------------------------
; Check on reproducibility of K04 figs

    l_sun = 3.826e33
    l_o2_raw = 10.D^s[sok].OII_3727_lum[0] * l_sun
    l_ha_raw = 10.D^s[sok].H_alpha_lum[0] * l_sun

    ebv_hahb = -2.5 * alog10(2.85/(s[sok].h_alpha[0]/s[sok].h_beta[0])) 
    ccm_orig_unred, [4861.33, 6562.82], [1.0, 1.0], 2.5, tenk
    k_ha = alog10(tenk[1])
    k_hb = alog10(tenk[0])
    s_e_bv = ebv_hahb / (k_hb[0] - k_ha[0]) > 0.02 ; following Kewley

    dered_fac = fltarr(4, nk)
    wl = [3727.0, 4861.3, 5006.8, 6562.8]
    for ii = 0, nk - 1 do begin
       ccm_orig_unred, wl, fltarr(4)+1, s_e_bv[ii], dered ; Pure CCM
       dered_fac[*,ii] = dered
    endfor
    l_ha = l_ha_raw * dered_fac[3,*]
    l_o2 = l_o2_raw * dered_fac[0,*]

; Note 4959 not included as per definition of K04!!!!
    r23 = alog10((s[sok].OII_3727[0] * dered_fac[0,*] + $
      s[sok].OIII_5007[0] * dered_fac[2,*]) / $
      (s[sok].H_beta[0] * dered_fac[1,*]))

; K04 equation 11 (ZKH94)
    oh_z94 = 9.265 - 0.33*r23 - 0.202*r23^2 - 0.207*r23^3 - 0.333*r23^4

;x = k.oh_z94
;x = s[sok].zstrong_12oh_zkh94
    x = oh_z94

    sfr_eq15 = 7.9e-42*l_o2 / (-1857.24 + 612.693*x - 67.0264*x^2 + 2.43209*x^3)
    sfr_ha = 7.9e-42 * l_ha


;----------------------
; plots

    set_plot, 'ps'
    device, filename = 'k04_test.ps', /inches

    empty = strarr(30) + ' '

    pagemaker, nx=2, ny=2, pos = pos, /normal, height = [2.0, 6.0], $
      ymargin = [1.0, 2.0], xmargin = [1.5, 0.5], xspace = [0.5]

; Kewley plot
    plot, k.sfr_ha_k98, alog10(k.sfr_ha_k98/k.sfr_o2_eq15), psym=4, /xlog, $
      xr=[0.005, 100], yr=[-0.5, 0.5], /xs, /ys, pos = pos[*,0], $
      xtickname = empty, title = 'Fig 12 b - K04', ytitle = 'log SFR ratio'
    oplot, [0.001, 100], [0, 0]

    plot, k.sfr_ha_k98, k.sfr_o2_eq15, psym=4, /xlog, /ylog, $
      xr = [0.005, 100], yr=[0.005, 100], /xs, /ys, $
      xtitle = 'SFR Ha k98', ytitle = 'SFR [OII] Egn 15', pos = pos[*,2], $
      /noerase
    oplot, [0.001, 100], [0.001, 100]
    xyouts, 0.01, 30, '!4r!X = ' + $
      string(robust_sigma(alog10(k.sfr_ha_k98/k.sfr_o2_eq15)), $
      format='(F5.3)'), charsize=1.5

; Now using ISPEC data
    plot, k.sfr_ha_k98, alog10(sfr_ha/sfr_eq15), psym=4, /xlog, $
      xr=[0.005, 100], yr=[-0.5, 0.5], /xs, /ys, pos = pos[*,1], /noerase, $
      xtickname = empty, ytickname = empty, title = 'Fig 12 b - ISPEC'
    oplot, [0.001, 100], [0, 0]

    plot, sfr_ha, sfr_eq15, psym=4, /xlog, /ylog, $
      xr = [0.005, 100], yr=[0.005, 100], /xs, /ys, $
      xtitle = 'SFR Ha k98', pos = pos[*, 3], ytickname = empty, $
      /noerase

    oplot, [0.001, 100], [0.001, 100]
    xyouts, 0.01, 30, '!4r!X = ' + $
      string(robust_sigma(alog10(sfr_ha/sfr_eq15)), format='(F5.3)'), $
      charsize=1.5

    device, /close
    set_plot, 'x'

return
end
