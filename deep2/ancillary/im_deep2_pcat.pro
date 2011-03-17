;+
; NAME:
;   deep2_pcat
; PURPOSE:
;   create comparison PCAT file to evaluate completeness
; CALLING SEQUENCE:
;   deep2_pcat
; REVISION HISTORY:
;   2005-04-07 MRB, NYU
;-
;------------------------------------------------------------------------------
pro im_deep2_pcat

pcatfiles=file_search(deep2_path(/analysis)+'pcat.*.fits.gz')
bigpcat=0
for i=0L, n_elements(pcatfiles)-1L do begin
    tmppcat=mrdfits(pcatfiles[i], 1)
    if(n_tags(bigpcat) eq 0) then $
      bigpcat=tmppcat $
    else $
      bigpcat=[bigpcat, tmppcat]
endfor

; split the photometric catalog into 45" chunks, and index to
; the first member of each group

pcat = bigpcat
ig=spheregroup(pcat.ra, pcat.dec, 1./3600., chunksize=45./3600., $
               firstgroup=firstgroup)
ii=where(firstgroup ne -1L, nii)
pcat=pcat[firstgroup[ii]]

; now cross-match against the spectroscopic catalog, and only keep
; photometric members within 30" of each spectroscopic target 

zcat=mrdfits(deep2_path(/analysis)+'zcat.dr3.uniq.good.fits.gz',1)

spherematch, zcat.ra, zcat.dec, pcat.ra, pcat.dec, 30./3600., $
  m1, m2, maxmatch=0

isort=sort(m2)
iuniq=uniq(m2[isort])
pcatm=pcat[m2[isort[iuniq]]]

mwrfits, pcatm, deep2_path(/analysis)+'deep2_pcat.fits', /create

return
end
