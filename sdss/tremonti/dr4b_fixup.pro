
dir = '/home/tremonti/sdss/samples/dr4_v5_1'

sfa = mrdfits( dir + '/gal_info_dr4_v5_1.fit', 1)
sfb = mrdfits( dir + 'b/gal_info_dr4_v5_1b.fit', 1)
nnew = n_elements(sfb)

ida = pmfspecid(s=sfa) 
idb = pmfspecid(s=sfb) 

indxa = cmset_op(ida, 'and', idb, /index)
indxb = cmset_op(idb, 'and', ida, /index)
;dr4_main = where(strmatch(sfb.release, 'dr*'))
sfa = 0
sfb = 0

junk = {junk: ' '}

files = [ $
         'gal_info_dr4_v5_1', $
         'gal_line_dr4_v5_1', $
         'gal_indx_dr4_v5_1', $
         'gal_misc_dr4_v5_1', $
         'gal_mass_dr4_v5_1', $
         'gal_oh_dr4_v5_1', $
         ;'gal_photo_dr4_fnal', $
         'photo_dr4_fnal_info', $
         'photo_dr4_fnal_match', $
         'photo_dr4_fnal_phot', $
         'photo_dr4_fnal_struct', $
         'sfr_fib', $
         'sfr_tot', $
         'specsfr_fib', $
         'specsfr_tot', $
         'tauV_res' $
         ]

for ii = 0, n_elements(files) -1  do begin
  print, files[ii]
  sa = mrdfits(dir + '/' + files[ii] + '.fit', 1)
  str = sa[0]
  struct_assign, junk, str  ; zero the elements

  sb = make_array(val = str, dim = nnew)
  sb[indxb] = sa[indxa]
  sa = 0
  mwrfits, sb, dir + 'b/' + files[ii] + '.fit', /create
  sb = 0
endfor

end
