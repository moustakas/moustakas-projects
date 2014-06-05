pro deep2_check_spec1d_dr4, debug=debug
; jm07aug30nyu - based on DEEP2_CHECK_SPEC1D.DR2 (but not backwards
;                compatible!); much, much faster and smarter
; jm08sep03nyu - some small tweaks; note that setting /ALL in
;   'read_deep2_zcat' results in NOINZCATINDX=-1; the default is to
;   read the 'uniq' redshift catalog
; jm13jun18siena - updated to DR4

    dr = 'dr4'

    datapath = deep2_path(/dr4)
    catalogs_path = deep2_path(/catalogs)

; read the full redshift catalog    
    zcat = mrdfits(catalogs_path+'zcat.deep2.dr4.uniq.fits.gz',1)
    ngalaxy = n_elements(zcat)

; first find all the 1D spectra on disk; replace some dumb EB, EC, and
; WB extensions to enable the matching
    spawn, 'find '+datapath+' -name "*.fits*" -print', spec1dfile_ondisk, /sh
    spec1dfile_ondisk_trim = repstr(repstr(repstr(spec1dfile_ondisk,'WB.',''),$
      'EB.',''),'EC.','')
    spec1dfile_ondisk_trim = file_basename(spec1dfile_ondisk_trim)
    
; now form the 1D spectrum name from the redshift catalog    
    slit = string(zcat.slit,format='(I3.3)')
    mask = string(zcat.mask,format='(I4.4)')
    obj = string(zcat.objname,format='(I0)')
    spec1dfile_zcat = 'spec1d.'+mask+'.'+slit+'.'+obj+'.fits.gz'
    
; remember, the CMSET_OP indices refer to the first argument; hence,
; ALLGOODINDX refers to SPEC1DFILE_ZCAT *not* SPEC1DFILE; you can
; check this by printing minmax(allgoodindx); also, all the CMSET_OP
; arrays need to be sorted
    notinzcatindx = cmset_op(spec1dfile_ondisk_trim,'and',/not2,spec1dfile_zcat,/index) ; spectrum exists, but not in ZCAT
    nospec1dindx = cmset_op(spec1dfile_zcat,'and',/not2,spec1dfile_ondisk_trim,/index)  ; in ZCAT, but no spectrum exists
    allgoodindx = cmset_op(spec1dfile_zcat,'and',spec1dfile_ondisk_trim,/index)         ; in ZCAT *and* spectrum exists
    allgoodindx_ondisk = cmset_op(spec1dfile_ondisk_trim,'and',spec1dfile_zcat,/index)  ; same as above, but for SPEC1DFILE_ONDISK

    spec1dfile_zcat_good = spec1dfile_zcat[allgoodindx]
    spec1dfile_ondisk_good = spec1dfile_ondisk[allgoodindx_ondisk] ; includes the full path name

; the only files that should fail this test is ones with EB, EC, and
; WB in the file name
;   verify = where(strtrim(spec1dfile_zcat_good,2) ne file_basename(spec1dfile_ondisk_good))
;   niceprint, spec1dfile_zcat_good[verify], file_basename(spec1dfile_ondisk_good[verify])
    help, spec1dfile_ondisk, spec1dfile_zcat, spec1dfile_zcat_good, $
      allgoodindx, allgoodindx_ondisk, nospec1dindx, notinzcatindx

    if (nospec1dindx[0] ne -1L) then begin
       splog, 'These should be equal: ', n_elements(spec1dfile_zcat), $
         n_elements(allgoodindx)+n_elements(nospec1dindx)
       zcat_nospec1d = zcat[nospec1dindx]
       splog, 'Writing '+catalogs_path+'zcat.'+dr+'.nospec1d.fits'
       mwrfits, zcat_nospec1d, catalogs_path+'zcat.'+dr+'.nospec1d.fits', /create
       spawn, 'gzip -f '+catalogs_path+'zcat.'+dr+'.nospec1d.fits', /sh
    endif

; look for duplicates - none in the 2007-Aug version of DR3!
;   for ii = 0L, ngalaxy-1L do begin
;      nobj = total(strmatch(spec1dfile,spec1dfile_zcat[ii]))
;      if (nobj gt 1.0) then splog, 'Multiple spectra for '+spec1dfile_zcat[ii]
;   endfor

;   good = mrdfits('../../analysis/zcat.dr3.uniq.good.fits.gz',1)
;   junk = mrdfits('../../analysis/zcat.dr3.junkspec1d.fits.gz',1)

    ngalaxy = n_elements(allgoodindx)
    moreinfo = replicate({galaxy: '', file: '', zcatindx: -1L, minwave: 0.0, $
      maxwave: 0.0, junkspec1d: 0, fluxflag: 0},ngalaxy)
;     maxwave: 0.0, junkspec1d: 0, extract_flag: 0},ngalaxy)
    zcat_out = struct_addtags(moreinfo,zcat[allgoodindx])
    zcat_out.file = strmid(spec1dfile_ondisk_good,strlen(datapath)) ; strip off DATAPATH (not FILE_BASENAME!)
    zcat_out.galaxy = repstr(file_basename(strtrim(zcat_out.file,2)),'.fits','')
    zcat_out.zcatindx = allgoodindx

; loop through to check for empty spectra and to store the minimum and
; maximum wavelengths
;   for ii = 34614L, ngalaxy-1L do begin
;   for ii = 15408, ngalaxy-1 do begin
    for ii = 0L, ngalaxy-1 do begin

       if (ii mod 500) eq 0 then print, 'Object ', ii
       print, format='("DEEP2_CHECK_SPEC1D: ",I0,"/",I0,".",A1,$)', $
         ii+1, ngalaxy, string(13b)
       print, ii

       spec1d = fill_gap(datapath+zcat_out[ii].file,/silent,header=hdr)
       
;      spec1 = mrdfits(datapath+zcat_out[ii].file,1,h1,/silent)
;      spec2 = mrdfits(datapath+zcat_out[ii].file,2,h2,/silent)

;; check to make sure that the blue and red extractions are right        
;       ex1 = strtrim(sxpar(h1,'EXTNAME'),2)
;       ex2 = strtrim(sxpar(h2,'EXTNAME'),2)
;       zcat_out[ii].extract_flag = (ex1 ne 'Bxspf-B') or (ex2 ne 'Bxspf-R')
       
;      zcat_out[ii].extname_blue = strmid(ex1,0,strpos(ex1,'-'))
;      zcat_out[ii].extname_red = strmid(ex2,0,strpos(ex2,'-'))
;      print, ex1, ex2
       
;       flux = [spec1.spec,spec2.spec]
;       ivar = [spec1.ivar,spec2.ivar]
;       wave = [spec1.lambda,spec2.lambda]

       if size(spec1d,/type) eq 4 or size(spec1d,/type) eq 2 then begin ; no data
          zcat_out[ii].junkspec1d = 1
       endif else begin
; make sure we can flux-calibrate
          cspec = qe_correction2(spec1d,hdr,pars,paramsendr,lambda,$
            throughtput_table,/telluric,ncorrect=7)
          spec = specphot(cspec,hdr,rlambda,rresponse,ilambda,iresponse)

          flux = spec.fluxl     ; [erg/s/cm2/A]
          ivar = spec.ivar_fluxl

;         flux = spec.spec
;         ivar = spec.ivar
          wave = spec.lambda
          fgood = where((ivar gt 0.0),nfgood)

          if total(finite(flux) eq 0) ne 0.0 or total(finite(ivar) eq 0) ne 0.0 then stop

;         if spec.fluxflag eq 3 then stop
          zcat_out[ii].fluxflag = spec.fluxflag
          if (nfgood eq 0L) then begin
;            zcat_out[ii].minwave = min(wave)
;            zcat_out[ii].maxwave = max(wave)
             zcat_out[ii].junkspec1d = 1
          endif else begin
             zcat_out[ii].minwave = min(wave[fgood])
             zcat_out[ii].maxwave = max(wave[fgood])

             if keyword_set(debug) then begin
                plot, wave, flux, ps=10, xsty=3, ysty=3 ;, xr=[7700,7850]
;               djs_oplot, spec1.lambda, spec1.spec, ps=10, color='blue'
;               djs_oplot, spec2.lambda, spec2.spec, ps=10, color='red'
                cc = get_kbrd(1)
             endif
          endelse
       endelse
    endfor 

;   ext = where(zcat_out.extract_flag,next)
;   if (next ne 0L) then for ie = 0L, next-1L do begin
;      h1 = headfits(datapath+zcat_out[ext[ie]].file,ext=1)
;      h2 = headfits(datapath+zcat_out[ext[ie]].file,ext=2)
;      print, zcat_out[ext[ie]].file, strtrim(sxpar(h1,'EXTNAME'),2), $
;        strtrim(sxpar(h2,'EXTNAME'),2)
;   endfor       

; write out    
    goodindx = where(zcat_out.junkspec1d eq 0 and zcat_out.fluxflag lt 3)
;   goodindx = where((zcat_out.junkspec1d eq 0) and (zcat_out.extract_flag eq 0))
    if (goodindx[0] ne -1L) then begin
       zcat_good = zcat_out[goodindx]
       splog, 'Writing '+catalogs_path+'zcat.'+dr+'.goodspec1d.fits'
       mwrfits, zcat_good, catalogs_path+'zcat.'+dr+'.goodspec1d.fits', /create
       spawn, 'gzip -f '+catalogs_path+'zcat.'+dr+'.goodspec1d.fits', /sh

; also write out the sample of Q>=3 objects
       q34 = where(zcat_good.zbest gt 0 and zcat_good.zquality ge 3)
       splog, 'Writing '+catalogs_path+'zcat.'+dr+'.goodspec1d.Q34.fits'
       mwrfits, zcat_good[q34], catalogs_path+'zcat.'+dr+'.goodspec1d.Q34.fits', /create
       spawn, 'gzip -f '+catalogs_path+'zcat.'+dr+'.goodspec1d.Q34.fits', /sh
    endif

; among the junk spectra, all of them have ZQUALITY<3; among the
; extract-flag objects, all but 8 have ZQUALITY<3
    junkindx = where(zcat_out.junkspec1d eq 1 or zcat_out.fluxflag ge 3)
;   junkindx = where(zcat_out.junkspec1d or zcat_out.extract_flag) 
    if (junkindx[0] ne -1L) then begin
       zcat_junk = zcat_out[junkindx]
       splog, 'Writing '+catalogs_path+'zcat.'+dr+'.junkspec1d.fits'
       mwrfits, zcat_junk, catalogs_path+'zcat.'+dr+'.junkspec1d.fits', /create
       spawn, 'gzip -f '+catalogs_path+'zcat.'+dr+'.junkspec1d.fits', /sh
    endif

stop    
    
return
end
    
pro deep2_check_junk

    datapath = deep2_path(/dr3)
    zcat = read_deep2_zcat(/good)
    ext = where(zcat.extract_flag,next)
    if (next ne 0L) then for ie = 0L, next-1L do begin
       h1 = headfits(datapath+strtrim(zcat[ext[ie]].file,2),ext=1)
       h2 = headfits(datapath+strtrim(zcat[ext[ie]].file,2),ext=2)
       print, zcat[ext[ie]].file, strtrim(sxpar(h1,'EXTNAME'),2)+', '+$
         strtrim(sxpar(h2,'EXTNAME'),2)
    endfor       

return    
end
