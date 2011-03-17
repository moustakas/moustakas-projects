;+
; NAME:
;   deep2_area
; PURPOSE:
;   get area of DEEP2 DR1
; CALLING SEQUENCE:
;   deep2_area, area
; REVISION HISTORY:
;   2005-11-01 MRB, NYU
;-
;------------------------------------------------------------------------------
pro im_deep2_area, area

pcatfiles=file_search(deep2_path(/analysis)+'pcat.*.fits.gz')
omegafield=0.48*0.72*(!DPI/180.)^2
omegatot=n_elements(pcatfiles)*omegafield
pcat=0
for i=0L, n_elements(pcatfiles)-1L do begin
    tmppcat=mrdfits(pcatfiles[i], 1,/silent)
    ii=where(tmppcat.pgal gt 0.2 and $
             deep2_color_cuts(tmppcat.magb,tmppcat.magr,tmppcat.magi) gt 0)
    tmppcat=tmppcat[ii]
    if(n_tags(pcat) eq 0) then $
      pcat=tmppcat $
    else $
      pcat=[pcat, tmppcat]
endfor
pcatm= mrdfits(deep2_path(/analysis)+'deep2_pcat.fits', 1,/silent)
ii=where(pcatm.pgal gt 0.2 and $
         deep2_color_cuts(pcatm.magb,pcatm.magr,pcatm.magi) gt 0)
pcatm=pcatm[ii]
area=float(n_elements(pcatm))/float(n_elements(pcat))*omegatot

end
