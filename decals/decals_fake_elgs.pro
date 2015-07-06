pro decals_fake_elgs, makestamps=makestamps, makeimages=makeimages, test=test
; jm14nov25siena - simulate 2D ELG images

    topdir = '/project/projectdirs/cosmo/work/decam/'
    fakedir = topdir+'cats/fake/'
    
    nfake = 100     ; number of fake galaxies to insert into the brick
    nstamp = 500L   ; number of fake ELGs to make
    nstampsize = 45 ; postage stamp width (symmetric) [pixels, roughly 12"]
    
    gain = 4.0       ; average DECam gain [e/ADU]
    pixscale = 0.262 ; average DECam pixel size [arcsec]
    
; choose the r-band magnitude distribution of the fake galaxies
    rbright = 18.0
    rfaint = 25.0 

; choose the bricks to muddle with fake galaxies
    bricks = '368723'
    
; --------------------------------------------------    
; build the images with the fake galaxies inserted
    if keyword_set(makeimages) then begin

; randomly choose NFAKE galaxies to build
       stampcat = mrdfits(fakedir+'stamps-cat.fits',1)
       nstamp = n_elements(stampcat)

       these = long(randomu(seed,nfake)*(nstamp-1))
       fake = stampcat[these]

; read the master zeropoints table, which we will need below
       zpts = mrdfits(topdir+'calib/photom/zeropoints.fits',1)
       
       brickinfo = mrdfits(topdir+'bricks.fits',1)
       for ii = 0, n_elements(bricks)-1 do begin

; get the ra,dec boundaries for this brick
          thisbrick = brickinfo[long(bricks)-1]
          ramin = thisbrick.ra1
          ramax = thisbrick.ra2
          decmin = thisbrick.dec1
          decmax = thisbrick.dec2

; randomly distribute the fake galaxies on the sky
          fake.ra = randomu(seed,nfake)*(ramax-ramin)+ramin
          fake.dec = randomu(seed,nfake)*(decmax-decmin)+decmin

; randomly assign r-band magnitudes and then gz magnitudes such that
; the grz colors are preserved
          fake.r = randomu(seed,nfake)*(rbright-rfaint)+rfaint
          fake.g = fake.r + fake.grz_true[0]-fake.grz_true[1]
          fake.z = fake.r + fake.grz_true[2]-fake.grz_true[1]
          
; get all the CCDs that touch this brick; to minimize overhead, loop
; on each unique CP image; it would be good to be able to generate
; this file on-the-fly
          ccdinfo = mrdfits(fakedir+'ccds-'+bricks[ii]+'.fits',1)
          allcpimage = strtrim(ccdinfo.cpimage,2)
          cpimage = allcpimage[uniq(allcpimage,sort(allcpimage))]
          
; generate all the directories we will need, write out the fake-galaxy
; catalog, and copy over the images which will be rewritten          
          cpdir = file_dirname(cpimage)
          cpdir = cpdir[uniq(cpdir,sort(cpdir))]
          for nn = 0, n_elements(cpdir)-1 do file_mkdir, fakedir+cpdir[nn]
          mwrfits, fake, fakedir+'fake-elgs.fits', /create
          write_ds9_regionfile, fake.ra, fake.dec, filename=fakedir+'fake-elgs.reg'
          
          ncpimage = n_elements(cpimage)
          for jj = 0, ncpimage-1 do begin
             thiscpimage = repstr(cpimage[jj],'.fz','')
             thiscpimage_ivar = repstr(thiscpimage,'ooi','oow') ; inverse variance map

             splog, 'Copying '+thiscpimage
             file_copy, topdir+'cats/'+thiscpimage, fakedir+thiscpimage, /overwrite
             file_copy, topdir+'cats/'+thiscpimage_ivar, fakedir+thiscpimage_ivar, /overwrite
             
             these = where(cpimage[jj] eq allcpimage,nccd)
             for kk = 0, nccd-1 do begin
                splog, 'Processing '+thiscpimage+', CCD'+strtrim(ccdinfo[these[kk]].cpimage_hdu,2)
                
                calname = strtrim(ccdinfo[these[kk]].calname,2)
                filt = strtrim(ccdinfo[these[kk]].filter,2)
                filttag = tag_indx(fake[0],filt)
             
; get the astrometry and the PSF for this CCD
                psffile = topdir+'calib/psfex/'+calname+'.fits'
                wcs = mrdfits(topdir+'calib/astrom/'+calname+'.wcs.fits',0,hdr,/silent)
                im = mrdfits(topdir+'cats/'+thiscpimage,ccdinfo[these[kk]].cpimage_hdu,imhdr,/silent)
                ivar = mrdfits(topdir+'cats/'+thiscpimage_ivar,ccdinfo[these[kk]].cpimage_hdu,/silent)
                ccdsz = size(im,/dim)
                
; get the magnitude zeropoint for this image
                match = where(file_basename(zpts.filename) eq file_basename(thiscpimage) and $
                  strtrim(zpts.ccdname,2) eq strtrim(ccdinfo[these[kk]].extname,2),nmatch)
                if nmatch ne 1 then message, 'Mismatch here!'
                thiszpt = zpts[match[0]]

; see common.py around line 372 for some info about the magnitude
; zeropoints 
                magzpt = thiszpt.ccdzpt + 2.5*alog10(thiszpt.exptime)
;               magzpt = sxpar(imhdr,'MAGZPT')

; just keep the objects that fall on this CCD; note that in order to
; use embed_stamp, the fake object has to be at least NSTAMPSIZE/2
; pixels from the edge (I think!)
                adxy, hdr, fake.ra, fake.dec, xx, yy
                keep = where(xx gt nstampsize/2.0 and xx lt ccdsz[0]-nstampsize/2.0-1 and $
                  yy gt nstampsize/2.0 and yy lt ccdsz[1]-nstampsize/2.0-1,ngal)
                
                if ngal gt 0 then begin
                   splog, 'Adding '+strtrim(ngal,2)+' fake galaxies'
                   galinfo = fake[keep]
                   galinfo.x = xx[keep]
                   galinfo.y = yy[keep]
                
; loop on each galaxy; read and normalize the stamp to the flux in
; this filter, convolve with the PSF, and then add noise
                   for gg = 0, ngal-1 do begin
                      splog, galinfo[gg].x, galinfo[gg].y, 'r=', galinfo[gg].(filttag)
                      stamp = mrdfits(fakedir+'stamps/stamp-'+$
                        string(galinfo[gg].stampid,format='(I4.4)')+'.fits',/silent)

                      stamp *= 10.0^(-0.4*galinfo[gg].(filttag))*1D9 ; -->nanomaggies
                      stamp *= 10.0^(0.4*(magzpt-22.5))              ; -->ADU
                      sz = size(stamp,/dim)

; blur by the PSF                      
                      psf = psfex_psf_recon(galinfo[gg].x,galinfo[gg].y,psffile)
                      blur_stamp = convolve(stamp,psf)

; add in Poisson noise; for simplicity assume a fixed gain of 4.0,
; which is the average of the values here:
; http://www.ctio.noao.edu/noao/sites/default/files/DECam/gain_rdnoise.txt
                      medivar = djs_median(ivar[(galinfo[gg].x-nstampsize/2.0)>0:$
                        (galinfo[gg].x+nstampsize/2.0)<(ccdsz[0]-1),$
                        (galinfo[gg].y-nstampsize/2.0)>0:$
                        (galinfo[gg].y+nstampsize/2.0)<(ccdsz[1]-1)])     ; [1/ADU^2]
                      blur_stamp_var = (blur_stamp*gain + gain^2/medivar) ; [electrons^2]

                      noisy_blur_stamp = (blur_stamp*gain + randomn(seed,nstampsize,nstampsize)*$
                        sqrt(blur_stamp_var))/gain ; [ADU]
;                     mwrfits, stamp, '~/junk.fits', /create
;                     mwrfits, blur_stamp, '~/junk.fits'
;                     mwrfits, noisy_blur_stamp, '~/junk.fits'

                      embed_stamp, im, noisy_blur_stamp, galinfo[gg].x-nstampsize/2.0, $
                        galinfo[gg].y-nstampsize/2.0               
                   endfor       ; close loop on GALAXY
                endif 
                djs_modfits, fakedir+thiscpimage, im, imhdr, exten_no=ccdinfo[these[kk]].cpimage_hdu

; -------------------------
; write out a couple test images                
                if keyword_set(test) and ngal gt 15 then begin
                   im0 = mrdfits(topdir+'cats/'+thiscpimage,0,hdr0)
                   mwrfits, im0, '~/junk.fits', hdr0, /create
                   mwrfits, im0, '~/junk2.fits', hdr0, /create

                   im = mrdfits(fakedir+thiscpimage,ccdinfo[these[kk]].cpimage_hdu,imhdr)
                   mwrfits, im, '~/junk.fits', imhdr

                   im = mrdfits(topdir+'cats/'+thiscpimage,ccdinfo[these[kk]].cpimage_hdu,imhdr)
                   mwrfits, im, '~/junk2.fits', imhdr
                   test = 0
                endif
; -------------------------
             endfor             ; close loop on CCDs
          endfor                ; close loop on CPIMAGE 
stop
       endfor                   ; close loop on BRICK
    endif
    
; --------------------------------------------------    
; build a fiducial set of galaxy images (stamps)     
    if keyword_set(makestamps) then begin
; read the parent sample of DEEP2 galaxies with morphological
; parameters; could potentially use a larger training set here from
; the ACS-GC; could also potentially restrict this to galaxies in the
; ELG/grz color box
       version = desi_elg_templates_version()
       templatepath = getenv('DESI_ROOT')+'/spectro/templates/'+$
         'elg_templates/'+version+'/'
       cat = mrdfits(templatepath+'elg_templates_obs_'+version+'.fits',1)

       keep = where(cat.radius_halflight gt 0 and cat.sersicn lt 6.0,ngal)
       cat = cat[keep]
       
; build a bootstrap mega-sample
       objindx = long(randomu(seed,nstamp)*(ngal-1))
    
       fake = {$
         stampid:    0L,$
;        filter:     '',$
         x:         0.0,$
         y:         0.0,$
         ra:         0D,$
         dec:        0D,$
         sersicn:   0.0,$
         r50:       0.0,$
         ba:        0.0,$
         phi0:      0.0,$
;        flux:      1.0,$
         g:         0.0,$
         r:         0.0,$
         z:         0.0,$
         grz_true: fltarr(3)}
       fake = replicate(fake,nstamp)
       fake.stampid = lindgen(nstamp)
       
;      fake.x = randomu(seed,nstamp)*(npix-1)
;      fake.y = randomu(seed,nstamp)*(npix-1)
       fake.sersicn = cat[objindx].sersicn
       fake.r50 = cat[objindx].radius_halflight/pixscale ; [pixel]
       fake.ba = cat[objindx].axis_ratio
       fake.phi0 = randomu(seed,nstamp)*360.0
       
; preserve the colors
       fake.grz_true[0] = cat[objindx].decam_g
       fake.grz_true[1] = cat[objindx].decam_r
       fake.grz_true[2] = cat[objindx].decam_z

       mwrfits, fake, fakedir+'stamps-cat.fits', /create
       
; build each stamp (normalized to unity) and then it write out
       for iobj = 0L, nstamp-1 do begin
          print, iobj, fake[iobj].ba, fake[iobj].sersicn, fake[iobj].r50, fake[iobj].phi0
          stamp = dfakegal(sersicn=fake[iobj].sersicn,r50=fake[iobj].r50,$
            flux=1.0,phi0=fake[iobj].phi0,ba=fake[iobj].ba,$
            nx=nstampsize,ny=nstampsize)
          mwrfits, stamp, fakedir+'stamps/stamp-'+string(iobj,format='(I4.4)')+'.fits', /create
       endfor
    endif
    
return
end
