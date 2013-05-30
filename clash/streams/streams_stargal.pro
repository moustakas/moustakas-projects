;+
; NAME:
;   dstargal_atlas
; PURPOSE:
;   identify stars and galaxies in NASA-Sloan Atlas image
; CALLING SEQUENCE:
;   dstargal_atlas 
; COMMENTS:
;   Requires dparents_atlas, dpsf_atlas to have been run.
;   Agressive about calling things stars (ie. tuned to large galaxies)
;   Outputs files in:
;     atlases/#/[base]-#-nimage.fits - star-subtracted images
;     atlases/#/[base]-#-sgset.fits - locations of stars, galaxies
; REVISION HISTORY:
;   11-Jan-2006  Written by Blanton, NYU
;-
;------------------------------------------------------------------------------
pro streams_stargal, plot=plot, gsmooth=gsmooth, glim=glim, gsaddle=gsaddle, $
  nsigma=nsigma, noclobber=noclobber

if(NOT keyword_set(ref)) then ref=1L
if(NOT keyword_set(gsmooth)) then gsmooth=10.
if(NOT keyword_set(glim)) then glim=25.
if(NOT keyword_set(gsaddle)) then gsaddle=50.
if(NOT keyword_set(nsigma)) then nsigma=20.
if(NOT keyword_set(maxnstar)) then maxnstar=3000L
if(NOT keyword_set(maxmem)) then maxmem=2.e+9

;; default to use base name same as directory name
spawn, 'pwd', cwd
base=(file_basename(cwd))[0]

;; read in pset
pset= gz_mrdfits(base+'-pset.fits',1)
imfiles=strtrim(pset.imfiles,2)
puse=pset.puse
nim= n_elements(imfiles)
nx=lonarr(nim)
ny=lonarr(nim)
images=ptrarr(nim)
nimages=ptrarr(nim)
ivars=ptrarr(nim)
hdrs=ptrarr(nim)
psfs=ptrarr(nim)

;; set up star and galaxy locations
sgset={base:base, $
       ref:ref, $
       iparent:-1L, $
       nstars:0L, $
       ra_stars:dblarr(maxnstar), $
       dec_stars:dblarr(maxnstar), $
       ngals:0L, $
       ra_gals:dblarr(maxnstar), $
       dec_gals:dblarr(maxnstar) }

;; find center object
phdr=gz_headfits(base+'-f160w.fits')
pim=gz_mrdfits(base+'-pimage.fits')
pcat=gz_mrdfits(base+'-pcat.fits',1)
npx=(size(pim,/dim))[0]
npy=(size(pim,/dim))[1]
iparent=pim[npx/2L, npy/2L]
xyad, phdr, float(npx/2L), float(npy/2L), raex, decex
cirrange, raex
sgset.iparent=iparent

if(iparent eq -1) then return

;; setup output directory
subdir= 'atlases'
spawn, /nosh, ['mkdir', '-p', subdir+'/'+strtrim(string(iparent),2)]

sgsetfile=subdir+'/'+strtrim(string(iparent),2)+'/'+base+'-'+ $
  strtrim(string(iparent),2)+'-sgset.fits'
if(gz_file_test(sgsetfile) ne 0 AND $
   keyword_set(noclobber) ne 0) then begin
    nimfile=subdir+'/'+strtrim(string(iparent),2)+'/'+base+ $
      '-nimage-'+strtrim(string(iparent),2)+'.fits'
    if(gz_file_test(nimfile) eq 0) then begin
        message, 'sgset file exists, but not nimfile'
    endif
    return
endif

;; read in the psf information
for k=0L, nim-1L do begin
   bimfile=(stregex(imfiles[k], '(.*)\.fits.*', /sub, /extr))[1]
   psfs[k]=ptr_new(mrdfits(bimfile+'-bpsf.fits'))
endfor

;; read in the images
ntest=10L
for k=0L, nim-1L do begin
   images[k]=ptr_new(gz_mrdfits('parents/'+base+'-parent-'+ $
                                strtrim(string(iparent),2)+'.fits',0+k*2L,hdr))
   hdrs[k]=ptr_new(hdr)
   nx[k]=(size(*images[k],/dim))[0]
   ny[k]=(size(*images[k],/dim))[1]
   ivars[k]=ptr_new(gz_mrdfits('parents/'+base+'-parent-'+ $
                               strtrim(string(iparent),2)+'.fits',1+k*2L))
endfor

;; find and subtract stars in the reference image
adxy, *hdrs[ref], raex, decex, xex, yex
nimages[ref]= ptr_new(streams_psfsub(image=*images[ref], ivar=*ivars[ref], $
  psf=*psfs[ref], x=xstar_r, y=ystar_r, $
  flux=fluxstar_r, nstars=nstars, exx=xex, $
  exy=yex, exr=1.5, nsigma=nsigma))
stop

if(nstars gt 0) then begin
   xyad, *hdrs[ref], xstar_r, ystar_r, ra_stars, dec_stars
   cirrange, ra_stars
   fluxstar= fltarr(n_elements(fluxstar_r), nim)
   fluxstar[*,ref]= fluxstar_r
endif

;; subtract same stars in all other images
if(nstars gt 0) then begin
   for k=0L, nim-1L do begin
      if(k ne ref) then begin
         adxy, *hdrs[k], ra_stars, dec_stars, xstar, ystar
         nimages[k]= ptr_new(dpsfsub_atlas(image=*images[k], ivar=*ivars[k], $
                                           psf=*psfs[k], x=xstar, y=ystar, $
                                           flux=tmp_fluxstar))
         fluxstar[*,k]=tmp_fluxstar
      endif
   endfor
endif else begin
   for k=0L, nim-1L do $
      if(k ne ref) then $
         nimages[k]= ptr_new(*images[k])
endelse

;; output subtracted images
nimfile=subdir+'/'+strtrim(string(iparent),2)+'/'+base+ $
  '-nimage-'+strtrim(string(iparent),2)+'.fits'
for k=0L, nim-1L do $
   mwrfits, *nimages[k], nimfile, *hdrs[k], create=(k eq 0)

;; find galaxies in reference image
simage=dsmooth(*nimages[ref], gsmooth)
ssig=dsigma(simage, sp=long(gsmooth*5.))
saddle=gsaddle*ssig
dpeaks, simage, xc=xc, yc=yc, sigma=ssig, minpeak=glim*ssig, npeaks=ngals, $
        /check, saddle= gsaddle, /refine

;; if no galaxies found, let's just assume there is one at center
if(ngals eq 0) then begin
    ngals=1L
    ra_gals=[raex]
    dec_gals=[decex]
    adxy, *hdrs[ref], ra_gals, dec_gals, xc, yc
endif

if(ngals gt 0) then begin
   xgals=float(xc)
   ygals=float(yc)
   
   ;; refine the centers
   drefine, *nimages[ref], xgals, ygals, smooth=2., $
            xr=xrgals, yr=yrgals, box=long(5)
   xyad, *hdrs[ref], xrgals, yrgals, ra_gals, dec_gals
   cirrange, ra_gals
endif

;; only store up to maxnstars stars
if(nstars gt maxnstar) then begin
    nstars=maxnstar
    ra_stars= ra_stars[0L:nstars-1L]
    dec_stars= dec_stars[0L:nstars-1L]
endif

;; store locations in sgset
sgset.nstars= nstars
if(nstars gt 0) then begin
   sgset.ra_stars[0:nstars-1]= ra_stars
   sgset.dec_stars[0:nstars-1]= dec_stars
endif
sgset.ngals= ngals
if(ngals gt 0) then begin
   sgset.ra_gals[0:ngals-1]= ra_gals
   sgset.dec_gals[0:ngals-1]= dec_gals
endif 

;; get memory use
maxpix= float(max(nx*ny))
memuse= 4.*maxpix*float(ngals)
if(memuse gt maxmem) then begin
    splog, 'LIMITING MEMORY USE'
    sgset.ngals= long(maxmem/(4.*maxpix))
    limitedmem=1
endif

;; output star and galaxy info
sgsetfile=subdir+'/'+strtrim(string(iparent),2)+'/'+base+'-'+ $
  strtrim(string(iparent),2)+'-sgset.fits'
dhdr= dimage_hdr()
if(keyword_set(limitedmem)) then $
  sxaddpar, dhdr, 'LIMMEM', 1, '1 if limited NGALS due to memory'
mwrfits, sgset, sgsetfile, dhdr, /create

if(keyword_set(plot)) then begin
   atv, *nimages[ref]
   if(ngals gt 0) then begin
      adxy, *hdrs[ref], sgset.ra_gals[0:ngals-1], sgset.dec_gals[0:ngals-1], $
            xg, yg
      atvplot, xg, yg, psym=4, color='blue', th=4
   endif
   if(nstars gt 0) then begin
      adxy, *hdrs[ref], sgset.ra_stars[0:nstars-1], sgset.dec_stars[0:nstars-1], $
            xs, ys
      atvplot, xs, ys, psym=4, color='red'
   endif
endif

heap_free, psfs
heap_free, images
heap_free, nimages
heap_free, ivars
heap_free, hdrs

return 
end
;------------------------------------------------------------------------------
