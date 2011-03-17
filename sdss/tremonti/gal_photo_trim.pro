
dir = '/d2/tremonti/sdss/samples/dr4_v5_1b/'

tot = 567486

ptags= ['RUN', 'CAMCOL', 'RERUN', 'FIELD', 'PARENT', 'ID', 'NCHILD', $
        'OBJC_TYPE', 'OBJC_PROB_PSF', 'OBJC_FLAGS', 'OBJC_FLAGS2', $
        'PSFCOUNTS', 'PSFCOUNTSERR', $
        'FIBERCOUNTS', 'FIBERCOUNTSERR', $
        'PETROCOUNTS', 'PETROCOUNTSERR', $ 
        'PETRORAD', 'PETRORADERR', $
        'PETROR50', 'PETROR50ERR', $
        'PETROR90', 'PETROR90ERR', $
        'ISO_A', 'ISO_AERR', 'ISO_B', 'ISO_BERR', $
        'R_DEV', 'R_DEVERR', 'AB_DEV', 'AB_DEVERR', $ 
        'R_EXP', 'R_EXPERR', 'AB_EXP', 'AB_EXPERR', $
        'COUNTS_MODEL', 'COUNTS_MODELERR', 'FRAC_PSF', $
        'RA', 'DEC', 'REDDENING', $
        'FIRSTMATCH', 'ROSATMATCH'] 


for ii = 0L, 283 do begin
  
  i1 = ii*2000
  if ii le 282 then i2 = (ii+1)*2000 - 1 else i2 = tot - 1
  print, ii, i1, i2

  pi = mrdfits(dir + 'gal_photo_dr4_fnal_v5_1c.fit', 1, range = [i1, i2])

  npi = struct_selecttags(pi,  select_tags=ptags)

  if ii eq 0 then np = npi else np = [np, npi]

endfor

mwrfits, np, dir + 'gal_photo_trim_dr4_v5_1c.fit', /create

end





