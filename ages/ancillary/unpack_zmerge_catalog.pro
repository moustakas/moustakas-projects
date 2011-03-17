; ###########################################################################
; RELEGATED - I AM NOW USING THE REDSHIFT FROM EISENSTEIN'S
;   KCORR CATALOG (see AGES_NDWFS_PHOTOMETRY) 
; ###########################################################################
;
;+
; NAME:
;       UNPACK_ZMERGE_CATALOG
;
; PURPOSE:
;       Parse the AGES/zmerge catalog.
;
; TODO:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2006 Feb 05, U of A - written
;-

pro unpack_zmerge_catalog, write=write, debug=debug
    
    catalogpath = ages_path(/catalogs)
    outpath = ages_path(/analysis)
    
    catalogfile = 'catalog.zmerge'
    fitsfile = catalogfile+'unpacked.fits'

; testing:

    targetcat = mrdfits(ages_path(/analysis)+'catalog.cat.noguidestars.fits.gz',1,/silent)
    if keyword_set(debug) then begin
       iband = mrdfits(ages_path(/analysis)+'catalog.ndwfsi.fits.gz',1,/silent)
       zmerge = mrdfits(ages_path(/analysis)+'catalog.zmerge.fits.gz',1,/silent)
       if (n_elements(zmerge) ne 0L) then goto, dosort
    endif
    
; initialize the output data structure

    zmerge0 = {$
      zmerge_catid:      0L,  $ ; row-ordered ID number
      zmerge_dmin:       0D0, $
      zmerge_ra:         0D0, $
      zmerge_dec:        0D0, $
      zmerge_ngood:      0L,  $
      zmerge_dgood:      0.0, $
      zmerge_snrtot:     0.0, $
      zmerge_z:          0.0, $
      zmerge_pass:       0L,  $
      zmerge_aper:       0L,  $
      zmerge_nduplicate: 0L,  $ ; number of duplicate entries in CATALOG.ZMERGE
      zmerge_nmultiple:  0L}    ; number of multiple (but not identical) entries in CATALOG.ZMERGE
    
; read an ASCII version of the full catalog and the sextractor header 

    splog, 'Reading '+catalogpath+catalogfile+'.'
    catalog = djs_readlines(catalogpath+catalogfile)
    ncat = n_elements(catalog)

    t0 = systime(1)
    for i = 0L, ncat-1L do begin

       print, format='("Object = ",I0,"/",I0,".",A4,$)', i+1L, ncat, string(13b)
       
       data = strsplit(catalog[i],' ',/extract) 
       ndata = n_elements(data)

       if (long(data[0]) eq 1L) and (long(data[5]) ge 1L) then begin ; positive match and at least one good redshift

          zmerge0.zmerge_catid  = i
          zmerge0.zmerge_dmin   = double(data[1])
          zmerge0.zmerge_ra     = double(data[2])
          zmerge0.zmerge_dec    = double(data[3])
          zmerge0.zmerge_ngood  = long(data[5])
          zmerge0.zmerge_dgood  = float(data[7])
          zmerge0.zmerge_snrtot = float(data[9])

          ngood = zmerge0.zmerge_ngood < 3L  ; only the first *three* good spectra are tabulated
          zmerge1 = replicate(zmerge0,ngood) ; replicate for each good redshift

          zindx = lindgen(ngood)*4L+10L
          pindx = lindgen(ngood)*4L+12L
          aindx = lindgen(ngood)*4L+13L

          if (max(zindx) ge ndata) or (max(pindx) ge ndata) or (max(aindx) ge ndata) then begin
             splog, 'Problem on line '+string(i,format='(I0)')+'.'
             stop
          endif
          
          zmerge1.zmerge_z    = float(data[zindx])
          zmerge1.zmerge_pass = long(data[pindx])
          zmerge1.zmerge_aper = long(data[aindx])

;         if (ngood ge 3L) then begin
;            struct_print, zmerge1
;            cc = get_kbrd(1)
;         endif

          if (n_elements(zmerge) eq 0L) then zmerge = zmerge1 else zmerge = struct_append(zmerge,zmerge1)
          
       endif

    endfor
    splog, 'Total time to unpack ZMERGE catalog = '+string(systime(1)-t0,format='(G0.0)')+' seconds.'

    dosort:

; now loop back through the catalog and for each unique combination of
; pass and aperture, throw away duplicates ZMERGE entries

    zmerge = zmerge[sort(zmerge.zmerge_pass)]
    nobject = n_elements(zmerge)

    finalindx = -1L
    skipindx = -1L
    
    counter = 0L

    for i = 0L, nobject-1L do begin

       match = where((zmerge[i].zmerge_pass eq zmerge.zmerge_pass) and $
         (zmerge[i].zmerge_aper eq zmerge.zmerge_aper),nmatch)
       skip = where(i eq skipindx,nskip)

       if (nmatch gt 1L) and (nskip eq 0L) then begin

          splog, 'Multiple entries ('+string(nmatch,format='(I0)')+') for PASS/APER '+$
            string(zmerge[i].zmerge_pass,format='(I0)')+'/'+string(zmerge[i].zmerge_aper,format='(I0)')+'.'
          zmerge[i].zmerge_nduplicate = nmatch

; if the duplicate entries are identical then just pick the zeroth
; entry, otherwise pick the object with the minimum ZMERGE_DMIN and
; flag these objects in ZMERGE

          if (total(zmerge[match].zmerge_dmin-zmerge[match[0]].zmerge_dmin) ne 0.0) then begin
             zmerge[i].zmerge_nmultiple = nmatch
             mindx = where(zmerge[match].zmerge_dmin eq min(zmerge[match].zmerge_dmin),nmindx,comp=rejindx)
             if (nmindx gt 1L) then message, 'Stop!!'
             finalindx = [finalindx,match[mindx]]
          endif else begin
             finalindx = [finalindx,match[0L]]
          endelse
;         help, finalindx
          skipindx = [skipindx,match]

          struct_print, zmerge[match]
          if keyword_set(debug) then begin
;            struct_print, targetcat[zmerge[match].zmerge_catid]
;            struct_print, struct_trimtags(iband[zmerge[match].zmerge_catid],$
;              select=['I_X_IMAGE','I_Y_IMAGE','I_FLAG_DUPLICATE','I_ALPHA_J2000','I_DELTA_J2000'])
;            stop
          endif

          if (n_elements(zmerge_dup) eq 0L) then zmerge_dup = zmerge[match] else $
            zmerge_dup = struct_append(zmerge_dup,zmerge[match])
          if (n_elements(targetcat_dup) eq 0L) then targetcat_dup = targetcat[zmerge[match].zmerge_catid] else $
            targetcat_dup = struct_append(targetcat_dup,targetcat[zmerge[match].zmerge_catid])
          counter = counter + 1L

       endif else if (nskip eq 0L) then begin
          finalindx = [finalindx,i]
       endif ;else stop

    endfor
;   print, counter

; index the final catalog

;   print, finalindx[1L:n_elements(finalindx)-1L]
    zmerge_final = zmerge[finalindx[1L:n_elements(finalindx)-1L]] ; offset from the -1L
    
    if keyword_set(write) then begin
       splog, 'Writing '+outpath+fitsfile+'.'
       mwrfits, zmerge_final, outpath+fitsfile, /create
       spawn, ['gzip -f '+outpath+fitsfile], /sh
    endif

    if keyword_set(write) then begin
       zmerge_dupfile = repstr(fitsfile,'.fits','.duplicates.fits')
       splog, 'Writing '+outpath+zmerge_dupfile+'.'
       mwrfits, zmerge_dup, outpath+zmerge_dupfile, /create
       spawn, ['gzip -f '+outpath+zmerge_dupfile], /sh

       targetcat_dupfile = repstr(fitsfile,'.fits','.duplicates.fits')
       splog, 'Writing '+outpath+targetcat_dupfile+'.'
       mwrfits, targetcat_dup, outpath+targetcat_dupfile, /create
       spawn, ['gzip -f '+outpath+targetcat_dupfile], /sh
    endif
       
stop    
    
return
end
    
