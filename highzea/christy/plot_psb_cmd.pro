
make_jpg = 0

if keyword_set(make_jpg) then begin
  dir = '/home/tremonti/sdss/samples/dr4_v5_1b/'
  sf = mrdfits(dir + 'gal_info_dr4_v5_1b.fit', 1)
  sm = mrdfits(dir + 'gal_misc_dr4_v5_1c.fit', 1)

  gr = sm.ugriz_absmag_kc01[1] - sm.ugriz_absmag_kc01[2]
  mi = sm.ugriz_absmag_kc01[3] + 5*alog10(0.7)  ; correct for h=1 in kcorr

  ok = where(sf.z lt 0.2 and sf.z gt 0.04)

  plot2dh, mi[ok], gr[ok], [-15,-25], [0.0, 1.2], 'cmd', psym=3, /xs, /ys, $
           xtitle = 'i!U0.1!N', ytitle = '(g - r)!U0.1!N', img = img, $
           xysz=[7.5, 6.0], nx = 250, ny=200
  stopprint
  tv, img
  write_jpeg, 'cmd.jpg', 230B - img
endif

read_jpeg, 'cmd.jpg', img

startprint, 'psb_cmd', xysz=[7.5, 6.0]

plotimage, img, imgxr=[-15,-25], imgyr=[0.0, 1.2], xr=[-17, -25], $
           xtitle = 'i!U0.1!N', ytitle = '(g - r)!U0.1!N', $
           charthick=2, charsize = 1.5

m = mrdfits('/home/tremonti/spectra/mmt/models/psb_models_bc03_z05_cmd.fit',1)

itau = 1
mod_i = reform(m.gri_01[2,*,itau]) + 10.5
mod_gr = reform(m.gri_01[0,*,itau] - m.gri_01[1,*,itau])
mod_age = m.burst_age/1e6

out_age = [25, 50, 100, 250, 500, 750, 1000, 1500, 3000, 5000]
linterp, mod_age, mod_i, out_age, out_i
linterp, mod_age, mod_gr, out_age, out_gr

oplot, mod_i, mod_gr, thick=2, color=!blue
oplot, out_i, out_gr, color=!blue, thick=2, psym=1

xyouts, out_i, out_gr, string(out_age, format='(I4)') + ' Myr', color=!blue, $
        charthick=2, charsize = 0.8

;---------------------------------------------
; Now add PSBs

; Read SDSS filters
fdir = '/home/tremonti/local/idlutils/data/filters/'
readcol, fdir + 'sdss_jun2001_g_atm.dat', gwl, gfl
readcol, fdir + 'sdss_jun2001_r_atm.dat', rwl, rfl
readcol, fdir + 'sdss_jun2001_i_atm.dat', iwl, ifl

gal = mrdfits(hizea_path(/analysis) + 'hizea_mmt_vfit.fit', 1)
b = mrdfits(hizea_path(/analysis) + 'hizea_burstfit.fit', 1)
ngal = n_elements(gal)
for ii = 0, ngal - 1 do begin
  hizea_join_data, gal[ii], restwl, rflux, rerr, mask
  redwl = (findgen(3500) + max(restwl) + 1)   
  bigwl = [restwl, redwl]
  cont = hizea_model_restore(m, b[ii], bigwl, gal[ii].vdisp)
;  bigfl = [rflux, cont[n_elements(rflux):*]]
;  bigfl[where(mask eq 1)] = cont[where(mask eq 1)]
  bigfl = cont

  ; only include galaxy light
  b[ii].plaw_ampl = 0
  b[ii].qso_ampl = 0
  gcont = hizea_model_restore(m, b[ii], bigwl, gal[ii].vdisp)

  gmag = filtermag(bigwl*1.1, bigfl*1d-17, gwl, gfl)
  rmag = filtermag(bigwl*1.1, bigfl*1d-17, rwl, rfl)
  imag = filtermag(bigwl*1.1, bigfl*1d-17, iwl, ifl)
  dist = lumdist(gal[ii].z)
  imag = imag - 5 * alog10(dist*1e6) + 5
  print, ii, gal[ii].short_name, imag

  oplot, [imag], [gmag - rmag], psym = 5, thick=3

  g_gmag = filtermag(bigwl*1.1, gcont*1d-17, gwl, gfl)
  g_rmag = filtermag(bigwl*1.1, gcont*1d-17, rwl, rfl)
  g_imag = filtermag(bigwl*1.1, gcont*1d-17, iwl, ifl)
  g_imag = g_imag - 5 * alog10(dist*1e6) + 5
  oplot, [imag, g_imag], [gmag - rmag, g_gmag-g_rmag], linestyle=0

  linterp, m.burst_age, indgen(24), b[ii].age + 1e9, aindx
  oaindx = where(m.burst_age eq b[ii].age)
  b[ii].age = m.burst_age[aindx]
  ocont = hizea_model_restore(m, b[ii], bigwl, gal[ii].vdisp)
  o_gmag = filtermag(bigwl*1.1, ocont*1d-17, gwl, gfl)
  o_rmag = filtermag(bigwl*1.1, ocont*1d-17, rwl, rfl)
  o_imag = filtermag(bigwl*1.1, ocont*1d-17, iwl, ifl)
  o_imag = o_imag - 5 * alog10(dist*1e6) + 5

  tindx = where(m.burst_tau eq b[ii].tau)
  di = m.gri_01[2,aindx,tindx] - m.gri_01[2,oaindx,tindx]

  oplot, [g_imag + di], [o_gmag - o_rmag], psym=5, thick=3, color=!red
  ;oplot, [g_imag, g_imag +di], [g_gmag - g_rmag, o_gmag-o_rmag], $
  ;        linestyle=0, color=!red

endfor

stopprint

end
