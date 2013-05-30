;+
; NAME:
;   dpsfselect_atlas
; PURPOSE:
;   check how well a PSF fits some data
; CALLING SEQUENCE:
;   ispsf= dpsfselect_atlas(image, ivar, x, y [, psf=, amp= ])
; INPUTS:
;   image - [nx, ny] input image
;   ivar - [nx, ny] input invverse variance
;   x, y - [N] positions to check
;   psf - [npx, npy] PSF image
; OPTIONAL INPUTS:
;   clipsize - size of region to calculation chi^2 from, in units
;              of the PSF FWHM (default 5)
; OUTPUTS:
;   ispsf - [N] 1 for PSF, 0 otherwise
;   amp - amplitude of fit (after peak-normalizing PSF)
; REVISION HISTORY:
;   1-Aug-2010  Written by Blanton, NYU
;-
;------------------------------------------------------------------------------
function streams_psfselect, image, ivar, x, y, amp=amp, psf=psf, flux=flux, $
                           dof=dof, subimage=subimage, noclip=noclip, $
                           clipsize=clipsize

if(NOT keyword_set(clipsize)) then clipsize=5.

nx=(size(image,/dim))[0]
ny=(size(image,/dim))[1]

npx=(size(psf,/dim))[0]
npy=(size(psf,/dim))[1]

subimage=image

dfit_mult_gauss, psf, 1, amp, psfsig, model=model, /quiet 
fwhm=psfsig*2.*sqrt(2.*alog(2.))

if(NOT keyword_set(noclip)) then begin
   ncx=(ceil(fwhm[0]*clipsize*0.5)*2L+1) < npx
   ncy=(ceil(fwhm[0]*clipsize*0.5)*2L+1) < npy
   xst= npx/2L-ncx/2L
   xnd= xst+ncx-1L
   yst= npy/2L-ncy/2L
   ynd= yst+ncy-1L
endif else begin
   xst=0L
   xnd=npx-1L
   yst=0L
   ynd=npy-1L
   ncx=npx
   ncy=npy
endelse

cc=fltarr(ncx,ncy)+1.
xx=reform(findgen(ncx)#replicate(1., ncy), ncx*ncy)/float(ncy)-0.5
yy=reform(replicate(1., ncx)#findgen(ncy), ncx*ncy)/float(ncx)-0.5
rr=sqrt((xx-ncx*0.5)^2+(yy-ncy*0.5)^2)

cenpsf= psf[xst:xnd, yst:ynd]

cmodel=fltarr(11,ncx*ncy)
cmodel[0,*]=reform(cenpsf/max(psf), ncx*ncy)
cmodel[1,*]=cc
cmodel[2,*]=xx
cmodel[3,*]=yy
cmodel[4,*]=xx^2
cmodel[5,*]=yy^2
cmodel[6,*]=xx*yy
cmodel[7,*]=xx^3
cmodel[8,*]=xx^2*yy
cmodel[9,*]=xx*yy^2
cmodel[10,*]=yy^3

amp=fltarr(n_elements(x))
flux=fltarr(n_elements(x))
chi2=fltarr(n_elements(x))
dof=fltarr(n_elements(x))+float(n_elements(cenpsf))-4.
for i=0L, n_elements(x)-1L do begin 
   
   ;; extract cutout
   cutout_image=fltarr(ncx,ncy) 
   cutout_ivar=fltarr(ncx,ncy) 
   embed_stamp, cutout_image, subimage, ncx/2L-x[i], ncy/2L-y[i] 
   embed_stamp, cutout_ivar, ivar, ncx/2L-x[i], ncy/2L-y[i] 
   cutout_ivar=cutout_ivar>0.
; check for too many masked pixels (usually on the edge)
   if total(cutout_ivar gt 0)/(ncx*ncy*1.0) lt 0.95 then begin 
      chi2[i] = 1E6
      continue 
   endif
   
   ;; refine position
   dprefine, cutout_image, cenpsf, ncx/2L, ncx/2L, xr=xr, yr=yr, $
             invvar=cutout_ivar
   x[i]=x[i]+xr-float(ncx/2L)
   y[i]=y[i]+yr-float(ncy/2L)
   
   ;; re-extract cutout
   cutout_image=fltarr(ncx,ncy) 
   cutout_ivar=fltarr(ncx,ncy) 
   embed_stamp, cutout_image, subimage, ncx/2L-x[i], ncy/2L-y[i] 
   embed_stamp, cutout_ivar, ivar, ncx/2L-x[i], ncy/2L-y[i] 
   cutout_ivar=cutout_ivar>0.

   hogg_iter_linfit, cmodel, reform(cutout_image, ncx*ncy), $
                     reform(cutout_ivar, ncx*ncy), coeffs, nsigma=10
   amp[i]=coeffs[0] 
   flux[i]=coeffs[0] 
   chi2[i]= total((reform(coeffs#cmodel, ncx*ncy)- $
                   reform(cutout_image, ncx*ncy))^2* $
                  reform(cutout_ivar, ncx*ncy))
   
   submodel= (reform(coeffs#cmodel, ncx,ncy)-reform(cutout_image, ncx,ncy))
;  if chi2[i] lt 100 then stop
   if(arg_present(subimage)) then begin
      subpsf= coeffs[0]*psf
      embed_stamp, subimage, -subpsf, x[i]-float(npx/2L), y[i]-float(npy/2L)
   endif

endfor

return, chi2

end
;------------------------------------------------------------------------------
