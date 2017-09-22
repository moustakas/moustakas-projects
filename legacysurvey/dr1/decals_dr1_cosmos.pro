pro decals_dr1_cosmos, matchbricks=matchbricks, matchcosmos=matchcosmos
; jm15jun09siena -

    dr1mdir = '/global/work/decam/release/dr1m/'

; preliminary stuff
    if keyword_set(matchbricks) then begin
       all = mrdfits('~/dr1/decals-bricks.fits',1)
       ff = file_search(dr1mdir+'tractor/???/tractor-*.fits')
       match, all.brickname, repstr(repstr(file_basename(ff),$
         'tractor-',''),'.fits',''), m1, m2
       dr1m = all[m1]
       mwrfits, dr1m, dr1mdir+'dr1m-bricks.fits', /create

       cosmos = mrdfits(getenv('IM_ARCHIVE_DIR')+$
         '/data/cosmos/cosmos_phot_20060103.fits.gz',1)
       spherematch, cosmos.ra, cosmos.dec, dr1m.ra, $
         dr1m.dec, 0.25, m1, m2
       mwrfits, dr1m[m2], dr1mdir+'dr1m-cosmos-bricks.fits', /create
    endif

    if keyword_set(matchcosmos) then begin
       dr1m = mrdfits(dr1mdir+'dr1m-cosmos-bricks.fits',1)
       for ii = 0, n_elements(dr1m)-1 do begin
          t1 = mrdfits(dr1mdir+'tractor/'+strmid(dr1m[ii].brickname,0,3)+$
            '/tractor-'+dr1m[ii].brickname+'.fits',1)
          if n_elements(tt) eq 0L then tt = t1 else tt = [tt,t1]
       endfor
       mwrfits, tt, dr1mdir+'dr1m-tractor.fits', /create

       cosmos = mrdfits(getenv('IM_ARCHIVE_DIR')+$
         '/data/cosmos/cosmos_phot_20060103.fits.gz',1)
       spherematch, cosmos.ra, cosmos.dec, tt.ra, tt.dec, 1D/3600, m1, m2
       mwrfits, tt[m2], dr1mdir+'dr1m-tractor-cosmos.fits', /create
       mwrfits, cosmos[m1], dr1mdir+'dr1m-cosmos.fits', /create
    endif

    cosmos = mrdfits(dr1mdir+'dr1m-cosmos.fits',1)
    tractor = mrdfits(dr1mdir+'dr1m-tractor-cosmos.fits',1)
    decals_to_maggies, tractor, mm, ii, /decam_grz
    grz = -2.5*alog10(mm)

    istar = where(cosmos.star and cosmos.blend_mask eq 0)
    tstar = where(strtrim(tractor[istar].type,2) eq 'PSF')
    texp = where(strtrim(tractor[istar].type,2) eq 'EXP')
    tdev = where(strtrim(tractor[istar].type,2) eq 'DEV')
    tcomp = where(strtrim(tractor[istar].type,2) eq 'COMP')
    help, istar, tstar, texp, tdev, tcomp

    good = where(finite(grz[2,*]))
    miss = where(strtrim(tractor[good].type,2) eq 'PSF' and $
      cosmos[good].star eq 0 and cosmos[good].blend_mask eq 0)
    im_plothist, grz[2,good], bin=0.1, xr=[15,25]
    im_plothist, grz[2,good[miss]], bin=0.1, /fill, /over
    
    plot, cosmos[good[miss]].z_mag, cosmos[good[miss]].z_mag-$
      grz[2,good[miss]], psym=3, xrange=[15,28]
    
    
    plot, cosmos[istar[tstar]].g_mag, cosmos[istar[tstar]].g_mag-$
      grz[0,istar[tstar]], psym=3, xrange=[15,28]
    djs_oplot, cosmos[istar[texp]].g_mag, cosmos[istar[texp]].g_mag-$
      grz[0,istar[texp]], psym=3, color='orange'
    djs_oplot, cosmos[istar[tdev]].g_mag, cosmos[istar[tdev]].g_mag-$
      grz[0,istar[tdev]], psym=3, color='blue'
    djs_oplot, cosmos[istar[tcomp]].g_mag, cosmos[istar[tcomp]].g_mag-$
      grz[0,istar[tcomp]], psym=3, color='green'
    
    
    
    
    
    
stop    
    
    
return
end
