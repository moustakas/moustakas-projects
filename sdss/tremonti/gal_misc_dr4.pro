 
dir = '/home/tremonti/sdss/samples/dr4_v5_1b/'
cspeed = 2.99792458d5

read = 1

if keyword_set(read) then begin
  print, 'Reading ...'
  sf = mrdfits(dir + 'gal_info_dr4_v5_1b.fit', 1)
  z = sf.z
  sf = 0
  m = mrdfits(dir + 'gal_mass_dr4_v5_1b.fit', 1)
  mass = m.mass
 
  p = mrdfits(dir + 'gal_photo_trim_dr4_v5_1c.fit', 1)  
 
  read = 0
endif

str = {zzmax: 0.0, dist: 0.0, infiber: 0.0, conc: 0.0, $
       r50z_kpc: 0.0, logmu: 0.0, uniq: 0, $
       kc00: fltarr(5), ugriz_absmag_kc00: fltarr(5), mb_mass: 0.0, $
       kc01: fltarr(5), ugriz_absmag_kc01: fltarr(5), $
       ubvri_absmag_vega_kc00: fltarr(5)}

sm = make_array(val=str, dim=n_elements(z))

m_r = p.petrocounts[2] - p.reddening[2]
m_abs = m_r - 5 * alog10(3e5 * z / 70.0) - 25
z_max = 70.0 / 3e5 * 10.D^((-m_abs + 17.77 - 25) / 5)
sm.zzmax = z / z_max

sm.dist = lumdist(z, H0=70, Omega_m=0.3, Lambda0=0.7)
sm.conc = p.petror90[2]/p.petroR50[2]
sm.infiber = 10.0^(-0.4*(p.fibercounts[2] - p.petrocounts[2]))

adist = sm.dist / (1 + z)^2 ; angular diameter dist
sm.r50z_kpc = adist*1.0e3 * p.petror50[4]/206265 ; physical half-light r in kpc
sm.logmu = mass + alog10(0.5 / (!PI * sm.r50z_kpc^2))

readcol, dir + 'parsed_duplicates150.txt', index, flag, format='(L,I)'
sm.uniq = flag

;-----------------------------
; Blantons kcorrect - v4!!

sm.kc00 = sdss_kcorrect(z, tsobj=p, band_shift=0.0, $
                        mass=mb_mass, absmag=absmag00)
sm.ugriz_absmag_kc00 = absmag00
sm.mb_mass = alog10(mb_mass)

sm.kc01 = sdss_kcorrect(z, tsobj=p, band_shift=0.1, absmag=absmag01)
sm.ugriz_absmag_kc01 = absmag01

kcbessel = sdss2bessell(z, tsobj=p, band_shift=0.0, absmag=bessell_absmag)
sm.ubvri_absmag_vega_kc00 = bessell_absmag

mwrfits, sm, dir + 'gal_misc_dr4_v5_1c.fit', /create


end
