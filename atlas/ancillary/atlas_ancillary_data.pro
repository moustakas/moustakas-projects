;+
; NAME:
;       ATLAS_ANCILLARY_DATA
;
; PURPOSE:
;       Compile ancillary data for the spectral atlas. 
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2002 May 14 - originally written
;       jm02-04 - many updates
;       jm05jul21uofa - updated again
;       jm08jan30nyu - updated to match ATLAS_NED_WEBGET
;-

pro atlas_ancillary_data, atlas, write=write

    version = atlas_version(/ancillary)
    outpath = atlas_path(/analysis)

    basicname = 'atlas_ned.fits.gz'
    photoname = 'atlas_ned_photo.fits.gz'
    distname = 'atlas_distances.fits.gz'
    diamname = 'atlas_diameters.fits.gz'

    leda = atlas_read_leda()
    
    write_ancillary_data, atlas, datapath=outpath, outpath=outpath, $
      basicname=basicname, photoname=photoname, distname=distname, $
      diamname=diamname, leda=leda

; set the redshift for this object to that of the redshift of the
; whole system (jm08apr14nyu)

    indx = speclinefit_locate(atlas,'ugc09425ned01')      ; =UGC09425
    refindx = speclinefit_locate(atlas,'ugc09425',/exact) ; =UGC09425NW
    atlas[indx].z = atlas[refindx].z
    
; for these objects the NED position is wrong

    indx = speclinefit_locate(atlas,'mrk0960') ; SIMBAD position; NED based on 2002PrivC.U..C....P
    atlas[indx].ra ='00:48:35.44'
    atlas[indx].dec = '-12:42:59.90'
    atlas[indx].position_ref = 'SIMBAD'

    indx = speclinefit_locate(atlas,'ic1623b') ; LEDA position; NED based on 2MASS
    atlas[indx].ra ='01:07:48.179'
    atlas[indx].dec = '-17:30:23.89'
    atlas[indx].position_ref = 'LEDA'

    indx = speclinefit_locate(atlas,'ugc05028') ; LEDA position; NED based on 2MASS
    atlas[indx].ra ='09:27:50.13'
    atlas[indx].dec = '+68:24:42.7'
    atlas[indx].position_ref = 'LEDA'

    indx = speclinefit_locate(atlas,'ugca225') ; LEDA position; NED based on 1986ApJS...61..305G
    atlas[indx].ra ='11:04:58.34'
    atlas[indx].dec = '+29:08:17.5'
    atlas[indx].position_ref = 'LEDA'

    indx = speclinefit_locate(atlas,'ugc06541') ; SIMBAD position; NED based on "NED measurement"
    atlas[indx].ra ='11:33:29.12'
    atlas[indx].dec = '+49:14:17.4'
    atlas[indx].position_ref = 'SIMBAD'

    indx = speclinefit_locate(atlas,'ngc5954') ; LEDA position; NED based on 2MASS
    atlas[indx].ra ='15:34:35.07'
    atlas[indx].dec = '+15:12:00.4'
    atlas[indx].position_ref = 'LEDA'

    indx = speclinefit_locate(atlas,'ngc5954') ; SIMBAD position; NED based on SDSS
    atlas[indx].ra = '23:17:36.61'
    atlas[indx].dec =  '+14:00:03.0'
    atlas[indx].position_ref = 'SIMBAD'
    
    indx = speclinefit_locate(atlas,'ngc7592a') ; LEDA position; no NED reference
    atlas[indx].ra = '23:18:21.73'
    atlas[indx].dec = '-04:24:57.4'
    atlas[indx].position_ref = 'LEDA'

; RC3 diameter is wrong ("NED essential note")

    indx = speclinefit_locate(atlas,'ngc0694') ; Hyperleda position; NED based on 2MASS
    atlas[indx].d25_maj =  0.1*10^0.76 ; [arcmin]
    atlas[indx].d25_min = atlas[indx].d25_maj*10.0^(-0.17) ; [arcmin]
    atlas[indx].d25_origin = 'LEDA'
    atlas[indx].d25_ref = 'LEDA'
    
; write out

    outname = 'atlas_ancillary_data_'+version+'.fits'
    
    if keyword_set(write) then begin
    
       splog, 'Writing '+outpath+outname+'.'
       mwrfits, atlas, outpath+outname, /create
       spawn, ['gzip -f '+outpath+outname], /sh

    endif

return
end

;;; compile literature metallicities by cross-matching with our HII
;;; region database
;;
;;    hii = read_hii_regions(hiined=bighiined)
;;    ra = 15.0*im_hms2dec(atlas.ra)
;;    dec = im_hms2dec(atlas.dec)
;;    
;;    hiira = 15.0D*im_hms2dec(bighiined.ra)
;;    hiidec = im_hms2dec(bighiined.dec)
;;
;;    ntot = im_djs_angle_match(ra,dec,hiira,hiidec,dtheta=5.0/3600.0,$
;;      units='degrees',mcount=mcount,mindx=mindx,mdist=mdist)
;;    match = where(mindx ne -1L,nmatch,comp=nomatch,ncomp=nnomatch)
;;    splog, 'Matched '+string(nmatch,format='(I0)')+' galaxies with HII-region data.'
;;
;;    hiined = bighiined[mindx[match]]
;;    
;;;   niceprint, atlas[match].galaxy, hiined.galaxy
;;    
;;    for imatch = 0L, nmatch-1L do begin
;;
;;       indx = where(strtrim(hiined[imatch].galaxy,2) eq strtrim(hii.hii_galaxy,2),nindx)
;;       if (nindx ne 0L) then begin
;;
;;; distinguish between objects with gradients (spirals) and dwarfs with
;;; just multiple HII-region measurements; define a "dwarf" as one that
;;; has no HII-region radii measured
;;
;;          if (min(hii[indx].radius) eq max(hii[indx].radius)) then begin ; dwarf!
;;
;;;            niceprint, hii[indx].hii_galaxy, hii[indx].hii_region, hii[indx].zt_log12oh, hii[indx].reference
;;
;;; compute the weighted average of all the measurements 
;;             
;;             good = where(hii[indx].zt_log12oh gt -900.0,ngood)
;;             if (ngood ne 0L) then begin
;;
;;                zt1 = hii[indx[good]].zt_log12oh
;;                zt1err = hii[indx[good]].zt_log12oh_err
;;                
;;                if (ngood eq 1L) then begin
;;                   zt = zt1
;;                   zterr = zt1err
;;                   ref = hii[indx[good]].reference
;;                endif else begin
;;                   zt    = total(zt1/zt1err^2)/total(1.0/zt1err^2)
;;                   zterr = 1.0/sqrt(total(1.0/zt1err^2))
;;                   ref = strjoin(strtrim(hii[indx[good]].reference,2),'; ')
;;                endelse
;;                     
;;                atlas[match[imatch]].lit_log12oh     = zt
;;                atlas[match[imatch]].lit_log12oh_err = zterr
;;                atlas[match[imatch]].lit_log12oh_ref = ref
;;
;;                niceprint, string(atlas[match[imatch]].galaxy,format='(A12)'), $
;;                  string(zt,format='(F12.4)'), string(zterr,format='(F12.4)'), $
;;                  string(ref,format='(A0)')
;;
;;             endif
;;             
;;          endif
;;
;;       endif 
;;       
;;    endfor 

