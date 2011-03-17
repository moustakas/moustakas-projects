
mdir = '/home/tremonti/models/UVBLUE/p05/'
fullname = file_search(mdir + '*flx')
name = file_basename(fullname)
steff = strmid(name, 1, 5)
slogg = strmid(name, 7, 2)
usteff = steff[uniq(steff, sort(steff))]

ifile = '/home/tremonti/models/padova/isochrones/iso_jc_z040s.dat'
readcol, ifile, iso_log_age, mass_i, mass_a, logl, iso_logte, iso_logg, $
         mbol, mu, mb, mv, mr, mi, mj, mh, mk, flum, $
         format = '(F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,D)'

help, flum

iso_te = 10.0^iso_logte
uiso_log_age = iso_log_age[uniq(iso_log_age)]
uiso_age = 10.0^uiso_log_age
nage = n_elements(uiso_age)

s = readUVBLUE(mdir + 't32000g40p05k2.flx')
npix = n_elements(s.lam)
pop_synth = {wave: s.lam, flux: fltarr(npix, nage), $
             age: uiso_age}

for aa = 0, nage - 1 do begin
  ok = where(iso_log_age eq uiso_log_age[aa], nok)

  flux = 0
  for ii = 1, nok - 1 do begin
    junk = min(abs(usteff - iso_te[ok[ii]]), tindx)
    str_te = string(usteff[tindx], format='(I5.5)')

    uslogg = slogg[where(steff eq usteff[tindx])]  
    junk = min(abs(uslogg - iso_logg[ok[ii]]), gindx)
    str_logg = string(uslogg[gindx], format='(I2.2)')  
  
    sfile = 't' + str_te + 'g' + str_logg + 'p05k2.flx'
    s = readUVBLUE(mdir + sfile)
 
    nstars = flum[ok[ii]] - flum[ok[ii-1]]
    print, 'Mass ', mass_i[ok[ii]], '  N_stars = ', nstars*1e6  
    nwl = where(s.lam gt 4400 and s.lam lt 4500)
    flux = flux + s.flam * 10.0^flum[ok[ii]]

    plot, s.lam, flux, xr=[2500, 4500], title = mass_i[ok[ii]]
  endfor
  pop_synth.flux[*,aa] = flux
endfor

mwrfits, pop_synth, '/home/tremonti/mmt/models/uvblue_synth.fit', /create

end
