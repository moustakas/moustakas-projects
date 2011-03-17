

mfile ='~/home/research/synthesis/bc03/bc03_padova1994_chab_z05_ssp.fit'
;mfile ='/home/tremonti/models/bc03/bc03_padova1994_chab_z05_ssp.fit'
ssp = mrdfits(mfile, 1)

b = im_read_bc03(isedpath='/home/tremonti/models/bc03/', $
                 isedfile = 'const_sfr_2msun_z050.ised', minwave=1106, $
                 maxwave=7000)

t = findgen(50) * 2e8  ; 10 Gyr
t[0] = 1e5
sfr = fltarr(50) + 2.0

;m = make_arb_sfh(t, sfr, use_wave=[1100, 7000], xr=[2500,4500], $
;                  modelfile = mfile)

linterp, b.age, findgen(220), t, aindx
aindx = round(aindx)

l_sun =  3.827d33 ; erg/s

for ii = 1, 9 do begin
   plot, b.wave, b.flux[*,aindx[ii]] * l_sun, title = b.age[aindx[ii]]/1e6, $
         xr=[3000, 7000], /xs

   cspec = add_past_sf(t, sfr, ssp, t[ii], wlrange=[1100,7000], wave=cwave)
   oplot, cwave, cspec * l_sun, color=!blue
 
   oplot, m.wave, m.flux[*,ii], color=!red
   print, b.age[aindx[ii]]
   wait, 2
endfor

end
