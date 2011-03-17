
pro mass_test, p=p, m=m

dir = '/home/tremonti/sdss/samples/mzpaper/'

if (not keyword_set(p) or not keyword_set(m)) then begin
  p = mrdfits(dir + 'mzgal_photo_dr2_mini_djs.fit', 1)
  m = mrdfits(dir + 'mzgal_mass_dr2.v5_0.fit', 1)
endif

g = p.petrocounts[1] - p.extinction[1]
r = p.petrocounts[2] - p.extinction[2]
gr= g - r 

g_kc = g - p.kcor_petro_00[1]
r_kc = r - p.kcor_petro_00[2]
gr_kc = g_kc - r_kc

nplothist, gr, bin=0.01, xr=[-1.5, 1.5]
nplothist, gr_kc, bin=0.01, /over, color=djs_icolor('red');!red

; From Blanton
B_ab = g_kc + 0.2747 + 0.4598 * (gr_kc - 0.7778)
B_vega = B_ab + 0.09

V_ab = r_kc + 0.2876 + 0.2246 * (gr_kc - 0.7778)
V_vega = V_ab -0.02

BV_vega = B_vega - V_vega


B_abs = B_vega + 5 - 5 * alog10(p.distance*1e6)
B_sun = 5.47 ; from BdJ 2001
L_B = 10.0^(-0.4 * (B_abs - B_sun))

; Bell & de Jong 2001
log_mlB_dietSP = -0.994 + 1.804 * BV_vega
log_mlB_kroupa = log_mlB_dietSP + 0.15
logm_BdJ = log_mlB_kroupa + alog10(L_B)


plot, m.mass, logm_BdJ, psym=3, xr=[6,12], yr=[6,12], $
      ytitle = 'log(M) Bell & de Jong', xtitle = 'log(M) Kauffmann'
oplot, [0, 20], [0, 20], color=!red

wait, 2

plothist, m.mass - logm_BdJ, bin=0.01, xr=[-1,1], $
  xtitle = 'Mass Offset (K - B)', ytitle = 'Number'

med = median(m.mass - logm_BdJ)
oplot, [0, 0], [0, 2000], color=!red
oplot, [1, 1] * med, [0, 2000], color=!green
print, 'Median Offset (K - B) ', med


end
