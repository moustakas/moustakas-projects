pro decals_fake_elgs, makestamps=makestamps
; jm14nov25siena - simulate 2D ELG images

    fakedir = getenv('HOME')+'/tmp/'

; read the parent sample of DEEP2 galaxies with morphological
; parameters 
    version = desi_elg_templates_version()
    templatepath = getenv('DESI_ROOT')+'/spectro/templates/'+$
      'elg_templates/'+version+'/'
    cat = mrdfits(templatepath+'elg_templates_obs_'+version+'.fits.gz',1)

    keep = where(cat.radius_halflight gt 0,ngal)
    cat = cat[keep]

; characteristics of each brick
    npix = 3600.0
    pixscale = 0.262 ; need to calculate this each time

    nx = 20 ; [pixels]
    ny = 20

    band = ['g','r','z']
    nband = n_elements(band)

; build a bootstrap mega-sample
    nfake = 2000L
    objindx = long(randomu(seed,nfake)*(ngal-1))

    rbright = 18.0
    rfaint = 25.0

    fake = {$
      x:         0.0,$
      y:         0.0,$
;     ra:         0D,$
;     dec:        0D,$
      sersicn:   0.0,$
      r50:       0.0,$
      ba:        0.0,$
      phi0:      0.0,$
      flux:      0.0,$
      grz:  fltarr(3)}
    fake = replicate(fake,nfake)

    fake.x = randomu(seed,nfake)*(npix-1)
    fake.y = randomu(seed,nfake)*(npix-1)
    fake.sersicn = cat[objindx].sersicn
    fake.r50 = cat[objindx].radius_halflight/pixscale ; [pixel]
    fake.ba = cat[objindx].axis_ratio
    fake.phi0 = randomu(seed,nfake)*360.0

; preserve the colors
    fake.grz[1] = randomu(seed,nfake)*(rfaint-rbright)+rfaint
    fake.grz[0] = fake.grz[1] + cat[objindx].decam_g-cat[objindx].decam_r
    fake.grz[2] = fake.grz[1] + cat[objindx].decam_z-cat[objindx].decam_r

    for iobj = 0L, nfake-1 do begin

; build the stamp and then normalize it to the appropriate flux
       stamp = dfakegal(sersicn=fake[iobj].sersicn,r50=fake[iobj].r50,$
         flux=1.0,phi0=fake[iobj].phi0,ba=fake[iobj].ba,$
         nx=nx,ny=ny)

       for ib = 1, 1 do begin
;      for ib = 0, nband-1 do begin
          stamp *= 10.0^(-0.4*fake[iobj].grz[ib]) ; [maggies]
          
; need the astrometry here and embed_stamp

       endfor       
         
stop
    endfor    

return
end
