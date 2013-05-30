;+
; NAME:
;   dfitpsf
; PURPOSE:
;   fit a PSF model to an image
; CALLING SEQUENCE:
;   dfitpsf, imfile
; INPUTS:
;   imfile - input FITS file
; OPTIONAL INPUTS:
;   natlas   - size of the PSF postage stamp (default 41)
;   maxnstar - maximum number of stars to use when constructing the
;              PSF (default 200)
;   
; KEYWORD PARAMETERS:
;   noclobber - if the PSF files already exist, do nothing and return 
; COMMENTS:
;   Currently always uses Nc=1 (no varying PSFs allowed)
;   Seems to work OK but many arbitrary parameters.
;   Ouputs (imfile is base.fits or base.fits.gz):
;     base-bpsf.fits - basic (single-fit) PSF
; REVISION HISTORY:
;   11-Jan-2006  Written by Blanton, NYU
;-
;------------------------------------------------------------------------------
pro streams_fitpsf, imfile, natlas=natlas, maxnstar=maxnstar, noclobber=noclobber, $
  cmap=cmap, base=base, seed=seed0, novpsf=novpsf, noivar=noivar, $
  check=check

if(NOT keyword_set(natlas)) then natlas=41L
if(NOT keyword_set(seed0)) then seed0=108L
seed=seed0

if(NOT keyword_set(base)) then $
  base=(stregex(imfile, '(.*)\.fits.*', /sub, /extr))[1]

if(keyword_set(noclobber)) then begin
    if(gz_file_test(base+'-bpsf.fits') gt 0 AND $
       gz_file_test(base+'-vpsf.fits') gt 0) then return
endif

splog, 'Reading image '+imfile
image=gz_mrdfits(imfile,/silent)
fits_info, imfile, n_ext=next, /silent
if (next ge 1L and keyword_set(noivar) eq 0) then begin
   splog, 'Reading inverse variance map'
   invvar = gz_mrdfits(imfile,1,/silent)
endif
nx=(size(image,/dim))[0]
ny=(size(image,/dim))[1]

;; set a bunch of parameters
plim=10.
box=natlas*2L
small=(natlas-1L)/2L
if (n_elements(maxnstar) eq 0L) then maxnstar=50L

; Set source object name
soname=filepath('libdimage.'+idlutils_so_ext(), $
                root_dir=getenv('DIMAGE_DIR'), subdirectory='lib')

;; median smooth the image and find and extract objects
xst=0L
yst=0L
msimage=dmedsmooth(image, invvar, box=box)
simage=image-msimage
dobjects, simage, objects=obj, plim=plim, seed=seed
dextract, simage, invvar, object=obj, extract=extract, small=small, $
  seed=seed
if(n_elements(extract) lt 3) then begin
   splog, 'not enough small enough objects in image'
   return
endif

extract1 = extract
splog, 'Identified ', n_elements(extract), ' candidate stars...'

; remove "stars" with more than 10% of the pixels masked (usually
; edges)
keep = (total(total(extract.atlas_ivar gt 0,2),1)/natlas^2.0) gt 0.95
keep = where(keep,nkeep)
if nkeep eq 0 then message, 'Rejected all stars!'
extract = extract[keep]
splog, '...keeping ', n_elements(extract), ' of them.'

rejsigma = [3.0,2.0,1.0,0.5]
for iter = 0L, n_elements(rejsigma)-1L do begin

   ;; do initial fit
   exatlas=fltarr(natlas*natlas, n_elements(extract))
   for i=0L, n_elements(extract)-1L do begin 
      scale=total(extract[i].atlas) 
      exatlas[*,i]=reform(extract[i].atlas/scale, natlas*natlas) 
   endfor
   bpsf=reform(djs_median( exatlas, 2), natlas, natlas)
   bpsf=bpsf/total(bpsf)

stop   
   
   ;; clip non-stars
   diff=fltarr(n_elements(extract))
   scale=fltarr(n_elements(extract))
   model=reform(bpsf, natlas*natlas)
   for i=0L, n_elements(extract)-1L do begin  
      scale[i]=total(model*reform(extract[i].atlas,natlas*natlas))/ $
               total(model*model) 
      diff[i]=total((extract[i].atlas-model*scale[i])^2*extract[i].atlas_ivar) 
   endfor
   
   ;; reject REJSIGMA outliers
   keep = where(diff lt mpchilim(rejsigma[iter],float(natlas)^2,/sigma),nkeep)
   splog, 'REJSIGMA = ', rejsigma[iter], '; keeping ', nkeep, ' stars'
   if (nkeep lt 3L) then begin
      splog, 'Not enough good stars found.'
      return
   endif
stop   
   extract = extract[keep]
   scale = scale[keep]
   
endfor

; find basic PSF
exatlas=fltarr(natlas*natlas, n_elements(extract))
model=reform(bpsf, natlas*natlas)
for i=0L, n_elements(extract)-1L do begin 
   scale[i]=total(model*reform(extract[i].atlas,natlas*natlas))/ $
            total(model*model) 
   exatlas[*,i]=reform(extract[i].atlas/scale[i], natlas*natlas) 
endfor
bpsf=reform(djs_median( exatlas, 2), natlas, natlas)
dfit_mult_gauss, bpsf, 1, amp, psfsig, model=model ; jm07may01nyu
bpsf=bpsf/total(model)

if(keyword_set(check)) then begin
   isort= reverse(sort(scale))
   for i=0L, n_elements(extract)-1L do begin
      mm= bpsf*total(model)*scale[isort[i]]
      atv, extract[isort[i]].atlas-mm, /block
   endfor
endif

pnx=(size(bpsf,/dim))[0]
pny=(size(bpsf,/dim))[1]
dfit_mult_gauss, bpsf, 1, amp, psfsig, model=model ; jm07may01nyu
mm=max(model)
gpsf=(model/mm) > 0.0001
if(keyword_set(novpsf)) then gpsf=fltarr(pnx,pny)+1.

; output basic PSF
mkhdr, hdr, 4, [natlas,natlas], /extend
sxdelpar, hdr, 'COMMENT' & sxdelpar, hdr, 'DATE'
sxaddpar, hdr, 'PSFSIGMA', float(psfsig[0]), ' Gaussian sigma [pixel]'
sxaddpar, hdr, 'NPSFSTAR', long(n_elements(extract)), $
          ' number of stars used in the PSF'
mwrfits, float(reform(bpsf, natlas, natlas)), base+'-bpsf.fits', hdr, /create
mwrfits, float(reform(model, natlas, natlas)), base+'-bpsf.fits', hdr 

npinit=3L
np=npinit
nc=3L
nchunk=2L
niter=3L

vatlas=fltarr(natlas*natlas, np^2*nchunk^2)
nb=lonarr(np^2*nchunk^2)
xx=fltarr(np^2*nchunk^2)
yy=fltarr(np^2*nchunk^2)
model=reform(bpsf, natlas*natlas)
for i=0L, np*nchunk-1L do begin
    for j=0L, np*nchunk-1L do begin
        ii=where(extract.xcen gt nx*(i)/(np*nchunk) AND $
                 extract.xcen le nx*(i+1)/(np*nchunk) AND $
                 extract.ycen gt ny*(j)/(np*nchunk) AND $
                 extract.ycen le ny*(j+1)/(np*nchunk), nii)
        if(nii gt 0) then begin
            sextract=fltarr(natlas*natlas, nii)
            for m=0L, nii-1L do begin
                sextract[*,m]= reform(extract[ii[m]].atlas, natlas*natlas)
                scale=total(model*sextract[*,m])/ total(model*model)
                sextract[*,m]=sextract[*,m]/scale
            endfor
            nb[i+(np*nchunk)*j]=nii	
            if(nii gt 1) then $
              vatlas[*,i+(np*nchunk)*j]= total(sextract, 2) $
            else $
              vatlas[*,i+(np*nchunk)*j]= sextract
            scale=total(model*vatlas[*,i+(np*nchunk)*j])/ total(model*model)
            vatlas[*,i+(np*nchunk)*j]= $
              (vatlas[*,i+(np*nchunk)*j]/scale - model)
            vatlas[*,i+(np*nchunk)*j]= vatlas[*,i+(np*nchunk)*j]*gpsf
            xx[i+(np*nchunk)*j]=mean(extract[ii].xcen)/float(nx)-0.5
            yy[i+(np*nchunk)*j]=mean(extract[ii].ycen)/float(ny)-0.5
        endif
    endfor
endfor

clipped=lonarr(np^2*nchunk^2)

splog, 'Iterating ...'
for iter=0L, niter-1L do begin
    iv=where(nb gt 0L and clipped eq 0, nv)
	  if(nv gt 4L) then begin
        if(nv gt npinit*npinit) then $
		      np=npinit $
        else $
          np=long(sqrt(nv))-1L

        em_pca, vatlas[*, iv], nc, eatlas, ecoeffs
        
        aa=dblarr(np*np, nv)
        k=0L
        for i=0L, np-1L do begin 
            for j=0, np-1L do begin 
                aa[k,*]=xx[iv]^(float(i))*yy[iv]^(float(j)) 
                k=k+1L 
            endfor 
        endfor
        
        xarr=(findgen(nx)#replicate(1.,ny))/float(nx)-0.5
        yarr=(replicate(1.,nx)#findgen(ny))/float(ny)-0.5
        
        coeffs=fltarr(np*np, nc)
        for c=0L, nc-1L do begin 
            sig=djsig(ecoeffs[c,*])
            weights=replicate(1., nv) /sig^2
            hogg_iter_linfit,aa, transpose(ecoeffs[c,*]), weights, $
              tmp_coeffs, nsigma=3., /true
            clipped[iv]=clipped[iv] OR (weights eq 0.)
            coeffs[*,c]=tmp_coeffs
        endfor

    endif else begin
        eatlas=fltarr(natlas, natlas, nc)
        coeffs=fltarr(np*np, nc)
    endelse
endfor


hdr=['']
sxaddpar, hdr, 'NP', np, 'number of polynomial terms'
sxaddpar, hdr, 'NC', nc, 'number of components in NMF'
sxaddpar, hdr, 'NATLAS', natlas, 'size of PSF image'
sxaddpar, hdr, 'XST', xst, 'start of source image used'
sxaddpar, hdr, 'YST', yst, 'start of source image used'
sxaddpar, hdr, 'NX', nx, 'dimension of source image (as used)'
sxaddpar, hdr, 'NY', ny, 'dimension of source image (as used)'
sxaddpar, hdr, 'SOFTBIAS', softbias, 'dimension of source image'
mwrfits, bpsf, base+'-vpsf.fits', hdr, /create
outatlas=reform(eatlas, natlas, natlas, nc)
for i=0L, nc-1L do $
  outatlas[*,*,i]=outatlas[*,*,i]/gpsf
mwrfits, outatlas, base+'-vpsf.fits'
outcoeffs=fltarr(200L, nc)
outcoeffs[0:np*np-1L, *]= coeffs
mwrfits, outcoeffs, base+'-vpsf.fits'

if(keyword_set(cmap)) then begin
    splog, 'Making cmap ...'
    out_cmap=fltarr(nx,ny,nc)
    if(min(coeffs) ne 0. OR max(coeffs) ne 0.) then begin
        for c=0L, nc-1L do begin
            k=0L 
            for i=0L, np-1L do begin 
                for j=0, np-1L do begin 
                    out_cmap[*,*,c]=out_cmap[*,*,c]+ $
                      coeffs[k,c]*xarr^(float(i))*yarr^(float(j)) 
                    k=k+1L 
                endfor 
            endfor 
        endfor
    endif
    mwrfits, out_cmap, base+'-vpsf.fits'
endif

end
;------------------------------------------------------------------------------
