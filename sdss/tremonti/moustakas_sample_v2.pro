
;******************************************************************************
; Read in data to make cuts
;******************************************************************************

dir = '/home/tremonti/sdss/samples/dr4_v5_1b/'

sf = mrdfits(dir + 'gal_info_dr4_v5_1b.fit', 1)
p1 = mrdfits(dir + 'photo_dr4_fnal_info_v5_1b.fit', 1)
p2 = mrdfits(dir + 'photo_dr4_fnal_photo_v5_1b.fit', 1)

;******************************************************************************
; Select subsets of data
;******************************************************************************

mainok = where(p2.petrocounts[2] - p1.reddening[2] gt 14.5 and $
               p2.petrocounts[2] - p1.reddening[2] le 17.77 and $
               (sf.sectarget and 8) eq 0 and $    ; not QA fiber
               (sf.primtarget and 64) ne 0 and $  ; targeted as main
               sf.z_warning eq 0) ; has a good redshift

infiber = 10.0^((p2.petrocounts[2] - p2.fibercounts[2])/2.5)

ok = cmset_op(mainok, 'and', where(sf.z gt 0.033 and infiber gt 0.10))

;---------------
; Get only DR2 plates 

sf_dr2 = mrdfits('/home/tremonti/sdss/samples/dr2_v5_0/' + $
                 'gal_info_dr2.v5_0.fit', 1)
dr2_plates = sf[uniq(sf_dr2.plateid)].plateid
sf_dr2 = 0

dr4_plates = sf[uniq(sf.plateid)].plateid

new_plates = cmset_op(dr4_plates, 'and', /not2, dr2_plates)

dr2 = ok
for ii = 0, n_elements(new_plates) - 1 do $
  dr2 = cmset_op(dr2, 'and', /not2, where(sf[dr2].plateid eq new_plates[ii]))

;clean duplicate plates

dr2 = cmset_op(dr2, 'and', where(strmatch(sf[dr2].release, 'dr1/2/3*')))

plothist, sf[ok].plateid
plothist, sf[dr2].plateid, /over, color=!red
help, dr2

;******************************************************************************
; Write out clipped sample
;******************************************************************************

outdir = '/home/tremonti/sdss/samples/moustakas/'

mwrfits, sf[dr2], outdir + 'sdss_info_dr2_v5_1b.fit', /create
sf = 0

mwrfits, p1[dr2], outdir + 'sdss_photo1_dr2_fnal.fit', /create
p1 = 0
mwrfits, p2[dr2], outdir + 'sdss_photo2_dr2_fnal.fit', /create
p2 = 0
p3 = mrdfits(dir + 'photo_dr4_fnal_struct.fit', 1)
mwrfits, p3[dr2], outdir + 'sdss_photo3_dr2_fnal.fit', /create
p3 = 0
p4 = mrdfits(dir + 'photo_dr4_fnal_match.fit', 1)
mwrfits, p4[dr2], outdir + 'sdss_photo4_dr2_fnal.fit', /create
p4 = 0

sl = mrdfits(dir + 'gal_line_dr4_v5_1b.fit', 1)
mwrfits, sl[dr2], outdir + 'sdss_line_dr2_v5_1b.fit', /create  
sl = 0

si = mrdfits(dir + 'gal_indx_dr4_v5_1b.fit', 1)
mwrfits, si[dr2], outdir + 'sdss_indx_dr2_v5_1b.fit', /create  
si = 0

m = mrdfits(dir + 'gal_mass_dr4_v5_1.fit', 1)
mwrfits, m[dr2], dir + 'sdss_mass_dr2_v5_1b.fit', /create
  
oh = mrdfits(dir + 'gal_oh_dr4_v5_1.fit', 1)
mwrfits, oh[dr2], dir + 'sdss_oh_dr2_v5_1b.fit', /create

sm = mrdfits(dir + 'gal_misc_dr4_v5_1b.fit', 1)
mwrfits, sm[dr2], outdir + 'sdss_misc_dr2_v5_1.fit', /create  
sm = 0

end

