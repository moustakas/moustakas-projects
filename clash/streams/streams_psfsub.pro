;+
; NAME:
;   dpsfsub_atlas
; PURPOSE:
;   subtract psfs 
; CALLING SEQUENCE:
;   subim= dpsfsub_atlas([imbase=, image=, ivar=, psf=, x=, y=, flux=])
; OPTIONAL INPUTS:
;   imbase - file base name (full name should be imbase+'.fits.gz')
;   image - [nx,ny] image to analyze
;   ivar - [nx,ny] ivar to analyze
;   psf - [npx,npy] PSF of image
;   nsigma - number of sigma to allow to call non-stellar (default 10)
;   xex, yex, rex - one coordinate and radius within which to exclude
;                   stars
; OPTIONAL INPUT/OUTPUTS:
;   x, y - [nstars] best fit location of stars
;   flux - [nstars] fluxes of each star
; OUTPUTS:
;   subim - attempt at fitting and subtracting PSFs
;   nstars - number of stars
; COMMENTS:
;   Either image, ivar, AND psf are set, or imbase is set
; REVISION HISTORY:
;   1-Mar-2006  Written by Blanton, NYU
;-
;------------------------------------------------------------------------------
function streams_psfsub, imbase=imbase, image=image, ivar=ivar, psf=psf, $
                        nsigma=nsigma, x=inout_x, y=inout_y, flux=flux, $
                        nstars=nstars, exx=xex, exy=yex, exr=rex

if(keyword_set(nsigma) eq 0) then $
   nsigma=20.

;; sanity checks and inputs
if(keyword_set(image) eq 0 AND $
   keyword_set(psf) eq 0 AND $
   keyword_set(ivar) eq 0) then begin
   if(keyword_set(imbase) eq 0) then $
      message, 'IMBASE must be set, or IMAGE, IVAR and PSF must be set'
   image=gz_mrdfits(imbase+'.fits',0)
   ivar=gz_mrdfits(imbase+'.fits',1)
   psf=mrdfits(imbase+'-bpsf.fits')
endif else begin
   if(keyword_set(image) eq 0 OR $
      keyword_set(psf) eq 0 OR $
      keyword_set(ivar) eq 0) then $
         message, 'IMBASE must be set, or IMAGE, IVAR and PSF must be set'
endelse

nx=(size(image,/dim))[0]
ny=(size(image,/dim))[1]

nix=(size(ivar,/dim))[0]
niy=(size(ivar,/dim))[1]
if(nix ne nx OR niy ne ny) then $
   message, 'IMAGE must be same size as IVAR'

;; find stars
nstars=0
if(keyword_set(inout_x) eq 0 OR keyword_set(inout_y) eq 0) then begin
   simplexy, image, xc, yc
   if(n_elements(xc) gt 0) then begin
      drefine, image, xc, yc, xr=x, yr=y
      chi2= streams_psfselect(image, ivar, x, y, psf=psf, dof=dof)

stop

      inex= bytarr(n_elements(x))
      if(n_elements(rex) gt 0) then begin
         inex= sqrt((x-xex)^2+(y-yex)^2) lt rex
      endif
      istar= where(chi2 lt dof+nsigma*sqrt(2.*dof) AND inex eq 0, nstars)
      if(nstars eq 0) then begin
         return, image
      endif 
      x=x[istar]
      y=y[istar]
      inout_x= x
      inout_y= y
   endif
endif else begin
   if(n_elements(inout_x) ne n_elements(inout_y)) then $
      message, 'X and Y must have same dimensions'
   if(n_elements(inout_x) gt 0) then begin
      xc= inout_x
      yc= inout_y
      nstars= n_elements(xc)
      drefine, image, xc, yc, xr=x, yr=y
   endif 
endelse 
if(nstars eq 0) then $
   return, image
flux=fltarr(nstars)

;; find approximate PSF size
dfit_mult_gauss, psf, 1, amp, psfsig, model=model, /quiet
fwhm=psfsig*2.*sqrt(2.*alog(2.))

;; set up model for subtraction
npx=(size(psf,/dim))[0]
npy=(size(psf,/dim))[1]

cc=fltarr(npx,npy)+1.
xx=reform(replicate(1., npx)#findgen(npy), npx*npy)/float(npx)-0.5
yy=reform(findgen(npx)#replicate(1., npy), npx*npy)/float(npy)-0.5
rr=sqrt((xx-npx*0.5)^2+(yy-npy*0.5)^2)

cmodel=fltarr(7,npx*npy)
cmodel[0,*]=reform(psf/total(psf), npx*npy)
cmodel[1,*]=cc
cmodel[2,*]=xx
cmodel[3,*]=yy
cmodel[4,*]=xx*xx
cmodel[5,*]=yy*yy
cmodel[6,*]=xx*yy

;; "subtract" each model
model=fltarr(nx,ny)
for i=0L, n_elements(x)-1L do begin 

   ;; cutout region around object center
   cutout_image=fltarr(npx,npy) 
   cutout_ivar=fltarr(npx,npy) 
   cutout_ivar=cutout_ivar>0. 
   embed_stamp, cutout_image, image, npx/2L-x[i], npy/2L-y[i] 
   embed_stamp, cutout_ivar, ivar, npx/2L-x[i], npy/2L-y[i] 
   
   ;; find best fit model including polynomial background
   hogg_iter_linfit, cmodel, reform(cutout_image, npx*npy), $
                     replicate(1., npx*npy), coeffs, nsigma=30
   flux[i]=coeffs[0] 
   tmodel= coeffs#cmodel
   
   ;; figure out what fraction of polynomial background the star is
   ;; (bounded between zero and one); subtract out that fraction of
   ;; the light of the full image
   fmodel= ((reform(coeffs[0]*cmodel[0,*]/tmodel,npx, npy))<1.)>0.
   pmodel= cutout_image*fmodel
   embed_stamp, model, pmodel, x[i]-float(npx/2L), y[i]-float(npy/2L)
   
endfor

nimage= image-model

return, nimage

end
;------------------------------------------------------------------------------
