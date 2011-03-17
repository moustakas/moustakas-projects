pro write_zintegrated_sample, intdust1, hii, write=write
; jm05sep09uofa - write integrated and nuc samples for the
;                 ZINTEGRATED paper
; jm06mar26uofa - updated

    sigmacut = 3.0
    
    outpath = atlas_path(/projects)+'zintegrated/'

    if (n_elements(intdust1) eq 0L) then intdust1 = read_integrated();atlasnodust=intnodust1)
    if (n_elements(hii) eq 0L) then hii = read_hii_regions()

; apply the S/N cuts on the integrated spectra

    keep = where($
;     (intdust1.nii_6584[0] /intdust1.nii_6584[1]  gt sigmacut) and $
;     (intdust1.oii_3727[0] /intdust1.oii_3727[1]  gt sigmacut) and $
;     (intdust1.oiii_5007[0]/intdust1.oiii_5007[1] gt sigmacut) and $
      (intdust1.h_alpha[0]  /intdust1.h_alpha[1]   gt sigmacut) and $
      (intdust1.h_beta[0]   /intdust1.h_beta[1]    gt sigmacut),nkeep)
    splog, 'Integrated spectra: '+string(nkeep,format='(I0)')+'/'+$
      string(n_elements(intdust1),format='(I0)')+' galaxies made the S/N cuts.'

    intdust = intdust1[keep]
;   intnodust = intnodust1[keep]
    
; cross-match against the HII-region database

    uhii = hii[uniq(hii.ned_galaxy,sort(hii.ned_galaxy))]
    match, strtrim(intdust.ned_galaxy,2), strtrim(uhii.ned_galaxy,2), indx1, indx2
;   niceprint, uhii[indx2].ned_galaxy, intdust[indx1].ned_galaxy

    match_uhii = uhii[indx2]
    doit = match_string(match_uhii.ned_galaxy,hii.ned_galaxy,index=index,/exact,/silent)
    match_hii = hii[index]
    nhii = n_elements(match_hii)
    
    match_intdust = intdust[indx1]
;   match_intnodust = intnodust[indx1]
    ngalaxy = n_elements(match_intdust)

    splog, 'Identified '+string(nhii,format='(I0)')+' HII regions in '+$
      string(ngalaxy,format='(I0)')+' galaxies.'
;   niceprint, match_uhii.ned_galaxy, match_intdust.ned_galaxy

; now figure out which of these galaxies have abundance gradients

    tempindx = lonarr(ngalaxy)-1L

    for i = 0L, ngalaxy-1L do begin

       indx = where(match_uhii[i].ned_galaxy eq match_hii.ned_galaxy,nindx)
       if (nindx lt 1L) then message, 'Hay problema.'

       if (min(match_hii[indx].hii_rc3_radius) ne max(match_hii[indx].hii_rc3_radius)) then begin
;         print, match_uhii[i].ned_galaxy, minmax(match_hii[indx].hii_rc3_radius)
          if (strmatch(match_uhii[i].hii_galaxy,'*NGC1068*') or strmatch(match_uhii[i].hii_galaxy,'*NGC1569*') or $
            strmatch(match_uhii[i].hii_galaxy,'*NGC4861*')) then begin
             if strmatch(match_uhii[i].hii_galaxy,'*NGC1068*') then splog, 'Rejecting NGC1068!' ; Seyfert 2
             if strmatch(match_uhii[i].hii_galaxy,'*NGC1569*') then splog, 'Rejecting NGC1569!' ; dwarf
             if strmatch(match_uhii[i].hii_galaxy,'*NGC4861*') then splog, 'Rejecting NGC4861!' ; dwarf
          endif else tempindx[i] = i
       endif
    endfor 
    
; subscript the final list, and cross-match with the parent HII-region
; database one last time
    
    finalindx = tempindx[where(tempindx ne -1L)]
    ngalaxy = n_elements(finalindx)

    final_intdust = match_intdust[finalindx]
;   final_intnodust = match_intnodust[finalindx]
    final_uhii = match_uhii[finalindx]
    
    doit = match_string(final_uhii.ned_galaxy,hii.ned_galaxy,index=index,/exact,/silent)
    final_hii = hii[index]
    nhii = n_elements(final_hii)

    splog, 'Retaining '+string(nhii,format='(I0)')+' HII regions in '+$
      string(ngalaxy,format='(I0)')+' galaxies with abundance gradients.'
;   srt = sort(final_intdust.galaxy)
;   niceprint, final_intdust[srt].galaxy, final_intdust[srt].ned_galaxy, final_uhii[srt].ned_galaxy

; NGC3351 and NGC4321 have S/N([O II]) < 3; set the flux equal to the
; upper limit and *artificially* set the uncertainty to the flux/3.0;
; in WRITE_ZINTEGRATED_TABLES and ZINTEGRATED we will distinguish
; these points as upper limits    

    losnr = where((final_intdust.oii_3727[0]/final_intdust.oii_3727[1] lt sigmacut),nlosnr)
    splog, 'Setting [O II] upper limits for '+strjoin(strtrim(final_intdust[losnr].galaxy,2),', ')+'.'
    final_intdust[losnr].oii_3727[1] = final_intdust[losnr].oii_3727[0]/3.01
    final_intdust[losnr].oii_3727_ew[1] = final_intdust[losnr].oii_3727_ew[0]/3.01

    final_intdust = struct_addtags(final_intdust,replicate({oh_lower_limit: 0L},ngalaxy))
    final_intdust[losnr].oh_lower_limit = 1L
    
; recompute the abundances and the extinction-corrected fluxes and
; abundances 

    zstrong = im_abundance(final_intdust,snrcut=sigmacut)
    final_intdust = struct_addtags(struct_trimtags(final_intdust,except=['ZSTRONG_*']),zstrong)

    final_intnodust = iunred_linedust(final_intdust,snrcut=sigmacut,/nopropagate)
    zstrong_nodust = im_abundance(final_intnodust,snrcut=sigmacut)
    final_intnodust = struct_addtags(struct_trimtags(final_intnodust,except=['ZSTRONG_*']),zstrong_nodust)
    
; now compute the dust ex    
    
    if keyword_set(write) then begin

       srt = sort(strtrim(final_intdust.galaxy,2))
       
       splog, 'Writing '+outpath+'zintegrated_intdust.fits.gz'
       mwrfits, final_intdust[srt], outpath+'zintegrated_intdust.fits', /create
       spawn, ['gzip -f '+outpath+'zintegrated_intdust.fits'], /sh

       splog, 'Writing '+outpath+'zintegrated_intnodust.fits.gz'
       mwrfits, final_intnodust[srt], outpath+'zintegrated_intnodust.fits', /create
       spawn, ['gzip -f '+outpath+'zintegrated_intnodust.fits'], /sh

       splog, 'Writing '+outpath+'zintegrated_hiiregions.fits.gz'
       mwrfits, final_hii, outpath+'zintegrated_hiiregions.fits', /create
       spawn, ['gzip -f '+outpath+'zintegrated_hiiregions.fits'], /sh

    endif

return
end
    
