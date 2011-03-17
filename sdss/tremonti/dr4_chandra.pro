
dir = '/home/tremonti/sdss/samples/dr4_v5_1b/'

st = {plate: 0, mjd: 0L, fiberid: 0, primtarget: 0L, $
     z: 0.0, ra: 0.0, dec: 0.0, distance: 0.0, infiber: 0.0, $
     OII_3727_cont: 0.0, OII_3727_cont_err: 0.0, $
     OII_3727_flux: 0.0, OII_3727_flux_err: 0.0, $    
     H_beta_cont: 0.0, H_beta_cont_err: 0.0, $    
     H_beta_flux: 0.0, H_beta_flux_err: 0.0, $    
     OIII_4959_cont: 0.0, OIII_4959_cont_err: 0.0, $    
     OIII_4959_flux: 0.0, OIII_4959_flux_err: 0.0, $    
     OIII_5007_cont: 0.0, OIII_5007_cont_err: 0.0, $    
     OIII_5007_flux: 0.0, OIII_5007_flux_err: 0.0, $    
     OI_6300_flux: 0.0, OI_6300_flux_err: 0.0, $    
     OI_6300_cont: 0.0, OI_6300_cont_err: 0.0, $    
     NII_6584_cont: 0.0, NII_6584_cont_err: 0.0, $
     NII_6584_flux: 0.0, NII_6584_flux_err: 0.0, $    
     H_alpha_cont: 0.0, H_alpha_cont_err: 0.0, $    
     H_alpha_flux: 0.0, H_alpha_flux_err: 0.0, $    
     SII_6717_cont: 0.0, SII_6717_cont_err: 0.0, $
     SII_6717_flux: 0.0, SII_6717_flux_err: 0.0, $    
     SII_6731_cont: 0.0, SII_6731_cont_err: 0.0, $
     SII_6731_flux: 0.0, SII_6731_flux_err: 0.0, $    
     extinction: fltarr(5), ba_exp: 0.0, $
     ugriz_petro: fltarr(5), ugriz_petro_err: fltarr(5), $
     ugriz_abs_petro: fltarr(5), ugriz_r50: fltarr(5), ugriz_r90:fltarr(5), $
     mass: 0.0, mass16: 0.0, mass84: 0.0, d4000_n_cor: 0.0, hda_cor: 0.0, $
     logmu: 0.0, release: '-', $
     log_OIII_lum: 0.0, class: 'unknown', d: 0.0, phi: 0.0, $
     v_disp: 0.0, v_disp_err: 0.0}

print, 'Reading Spectro Info'
sf = mrdfits(dir + 'gal_info_dr4_v5_1b.fit', 1)
ng = n_elements(sf)
s = make_array(val=st, dim=ng)

struct_assign, sf, s, /nozero
s.plate = sf.plateid
z = sf.z ; save for later
sf = 0
s.distance = lumdist(s.z, H0=70, Omega_m=0.3, Lambda0=0.7)

print, 'Reading Line Info'
sl = mrdfits(dir + 'gal_line_dr4_v5_1b.fit', 1)
struct_assign, sl, s, /nozero

; Combine OII lines
ok = where(sl.OII_3726_flux_err gt 0 and sl.OII_3729_flux_err gt 0)
s[ok].OII_3727_flux = sl[ok].OII_3726_flux + sl[ok].OII_3729_flux
s[ok].OII_3727_flux_err = sqrt(sl[ok].OII_3726_flux_err^2 + $
                               sl[ok].OII_3729_flux_err^2)
bad = where(sl.OII_3726_flux_err le 0 or sl.OII_3729_flux_err le 0)
s[bad].OII_3727_flux_err = -1
s.OII_3727_cont = (sl.OII_3726_cont + sl.OII_3729_cont) / 2
s.OII_3727_cont_err = sqrt(sl.OII_3726_cont_err^2 + sl.OII_3729_cont_err^2) / 2
sl = 0

print, 'Reading Masses'
m = mrdfits(dir + 'gal_mass_dr4_v5_1b.fit', 1)
struct_assign, m, s, /nozero
m = 0
; Re-do this so bad z's from mass file get overwritten.
s.z = z

; add Jarle masses
m = mrdfits(dir + 'sfr_info.fit', 1)
ok = where(s.mass ne 0, comp=nogmass)
plot, s[ok].mass, m[ok].lgm, psym=3
s[nogmass].mass = m[nogmass].lgm
m = 0

; These overwrite values in mass file (not available for all gal)
si = mrdfits(dir + 'gal_indx_dr4_v5_1b.fit', 1)
s.d4000_n_cor = si.d4000_n_sub
s.hda_cor = si.lick_hd_a_sub   
si = 0

print, 'Reading Photometry'
p = mrdfits(dir + 'photo_dr4_fnal_info_v5_1b.fit', 1)
struct_assign, p, s, /nozero
s.extinction = p.reddening
p = 0

p = mrdfits(dir + 'photo_dr4_fnal_struct_v5_1b.fit', 1)
s.ba_exp = p.ab_exp[2]
p = 0

p = mrdfits(dir + 'photo_dr4_fnal_photo_v5_1b.fit', 1)
struct_assign, p, s, /nozero
s.ugriz_petro = p.petrocounts
s.ugriz_petro_err = p.petrocountserr
s.infiber = 10.0^(-0.4*(p.fibercounts[2] - p.petrocounts[2]))
s.ugriz_r50 = p.petror50
s.ugriz_r90 = p.petror90
p = 0

dm = -5*alog10(s.distance*1e6) + 5
s.ugriz_abs_petro = s.ugriz_petro - s.extinction + ((fltarr(5)+1) # dm) 

adist = s.distance / (1 + s.z)^2 ; angular diameter dist
r50z_kpc = adist*1.0e3 * s.ugriz_r50[4]/206265 ; physical half-light r in kpc
s.logmu = s.mass + alog10(0.5 / (!PI * r50z_kpc^2))

; Compute AGN Params
s.log_OIII_lum = alog10(s.OIII_5007_flux * 1e-17 * $
                        4*!PI * (s.distance*1d6 * 3.086d18)^2)

sey = agnflag(s, sig = 7.0, n2ha = n2ha, o3hb = o3hb, gal = gal)
s[sey].class = 'seyfert2'
s[gal].class = 'galaxy'
x = n2ha[sey] + 0.45
y = o3hb[sey] + 0.5
s[sey].d = sqrt(x^2 + y^2)
s[sey].phi = atan(x/y) * !RADEG


;Don't include duplicates
readcol, dir + 'parsed_duplicates150.txt', index, dup, format='(L,I)'
ok = where(dup eq 1)
help, s, index, ok

mwrfits, s[ok], dir + 'gal_chandra_dr4_v5_1b.fit', /create

end
 
