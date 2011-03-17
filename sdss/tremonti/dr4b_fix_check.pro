

if keyword_set(read) then begin

  dir = '/home/tremonti/sdss/samples/dr4_v5_1b/'
  sla = mrdfits(dir + 'gal_line_dr4_v5_1.fit', 1)
  slb = mrdfits(dir + 'gal_line_dr4_v5_1b.fit', 1)

  sfb = mrdfits( dir + 'gal_info_dr4_v5_1b.fit', 1)
  sfa = mrdfits( dir + 'gal_info_dr4_v5_1.fit', 1)

  read = 0
endif

a = where(sfa.plateid ne 0)

plot, sfb[a].plateid, sfa[a].z/sfb[a].z, psym =3, yr=[0.5, 1.5]

wait, 2

a = where(sla.OIII_5007_flux_err > 0 and $
          sla.OIII_5007_flux/sla.OIII_5007_flux_err gt 7)

a4 = cmset_op(a, 'and', where(strmatch(sfb.release, 'dr*4*')))

o3cr = sla.OIII_5007_cont / slb.OIII_5007_cont
o3fr = sla.OIII_5007_flux / slb.OIII_5007_flux
o3sr = sla.sigam_forbidden / slb.sigma_forbidden

plot, slb[a].plateid, o3cr[a], psym=3, yr=[0, 10]

end
