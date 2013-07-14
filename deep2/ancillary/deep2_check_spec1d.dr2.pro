pro deep2_check_spec1d
; jm07feb03nyu

    datapath = deep2_path(/dr2)
    analysis_path = deep2_path(/analysis)

    zcat = read_deep2_zcat()
    ngalaxy = n_elements(zcat)

    spawn, 'find '+datapath+' -name "*.fits" -print', spec1dfile, /sh

;   for ii = 26670L, ngalaxy-1L do begin
    for ii = 0L, ngalaxy-1L do begin

       print, format='("DEEP2_CHECK_SPEC1D: ",I0,"/",I0,".",A1,$)', ii+1, ngalaxy, string(13b)

       slit = string(zcat[ii].slitname,format='(I3.3)')
       obj = string(zcat[ii].objno,format='(I0)')

       match = where(strmatch(spec1dfile,'*spec1d.????.*'+slit+'.'+obj+'.fits') eq 1B,nobj)
       if total(strmatch(spec1dfile[match],'*2002dec30/2002oct03*') ge 1.0) then begin ; ignore the *duplicate* data in the 2002dec30 subdirectory
          keep = where(strmatch(spec1dfile[match],'*2002dec30/2002oct03*') eq 0B,nobj)
          match = match[keep]
       endif

       if (nobj gt 1L) then begin ; grab the right one
          date = strsplit(strmid(date_conv(zcat[ii].date,'S'),0,11),'-',/extract)
          date = date[2]+strlowcase(date[1])+string(date[0],format='(I2.2)')
          thismatch = where(strmatch(spec1dfile[match],'*'+date+'*'),nobj)
          if (nobj ne 1L) then begin ; pick the object with the EB, EC, or WB in the filename (arbitrary!?!) 
             this = where(strmatch(spec1dfile[match[thismatch]],'*EB*') or strmatch(spec1dfile[match[thismatch]],'*EC*') or $
               strmatch(spec1dfile[match[thismatch]],'*WB*'),nobj)
             if (nobj ne 1L) then message, 'This should not happen!' else match = match[thismatch[this]]
          endif else match = match[thismatch]
       endif
       
       case nobj of
          0L: begin
             splog, 'No spectrum found for object ', ii
             if (n_elements(notfoundindx) eq 0L) then notfoundindx = ii else notfoundindx = [notfoundindx,ii]
          end
          1L: begin
             thisfile = strmid(spec1dfile[match[0]],strlen(datapath)) ; strip off the root path
             spec1 = mrdfits(spec1dfile[match[0]],1,/silent)
             spec2 = mrdfits(spec1dfile[match[0]],2,/silent)
             flux = [spec1.spec,spec2.spec]
             ivar = [spec1.ivar,spec2.ivar]
             wave = [spec1.lambda,spec2.lambda]
             fgood = where(ivar gt 0.0,nfgood)
             if (nfgood ne 0L) then begin
                if (n_elements(goodindx) eq 0L) then goodindx = ii else goodindx = [goodindx,ii]
                moreinfo1 = {file: thisfile, zcatindx: ii, minwave: min(wave[fgood]), maxwave: max(wave[fgood])}
                if (n_elements(moreinfo) eq 0L) then moreinfo = moreinfo1 else moreinfo = [moreinfo,moreinfo1]
             endif else begin 
                if (n_elements(junkindx) eq 0L) then junkindx = ii else junkindx = [junkindx,ii] ; these objects have no data, but also ZQUALITY<3
                if (n_elements(junkinfo) eq 0L) then junkinfo = {file: thisfile, zcatindx: ii} else $
                  junkinfo = [junkinfo,{file: thisfile, zcatindx: ii}]
             endelse
          end 
          else: begin
             message, 'Multiple matches for object ', ii, nobj ; this should never happen
          end
       endcase
       
    endfor

    if (goodindx[0] ne -1L) then begin
       zcat_good = struct_addtags(zcat[goodindx],moreinfo)
       splog, 'Writing '+analysis_path+'zcat.dr2.uniq.good.fits'
       mwrfits, zcat_good, analysis_path+'zcat.dr2.uniq.good.fits', /create
       spawn, 'gzip -f '+analysis_path+'zcat.dr2.uniq.good.fits', /sh
    endif

    if (junkindx[0] ne -1L) then begin
       zcat_junk = struct_addtags(zcat[junkindx],junkinfo)
       splog, 'Writing '+analysis_path+'zcat.dr2.junk.fits'
       mwrfits, zcat_junk, analysis_path+'zcat.dr2.junk.fits', /create
       spawn, 'gzip -f '+analysis_path+'zcat.dr2.junk.fits', /sh
    endif

;   if (notfoundindx[0] ne -1L) then begin
;      zcat_notfound = zcat[notfoundindx]
;      splog, 'Writing '+analysis_path+'zcat.dr2.notfound.fits'
;      mwrfits, zcat_notfound, analysis_path+'zcat.dr2.notfound.fits', /create
;      spawn, 'gzip -f '+analysis_path+'zcat.dr2.notfound.fits', /sh
;   endif
    
return
end
    
