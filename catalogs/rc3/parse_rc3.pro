pro parse_rc3, silent=silent
; jm03may4uofa

    rc3_path = filepath('',root_dir=getenv('CATALOGS_DIR'),subdirectory='rc3')

    pushd, rc3_path

    rc3file = findfile('rc3.fits*',count=fcount)
    if not keyword_set(silent) then splog, 'Reading the original RC3 galaxy catalog.'
    case fcount of
       0L: begin
          print, 'RC3 FITS file not found.'
          return
       end
       1L: begin
          rc3file = rc3file[0]
          if strmatch(rc3file,'*.gz') eq 1B then begin
             spawn, ['gunzip '+rc3file], /sh
             rc3fits = strmid(rc3file,0,strpos(rc3file,'.gz'))
          endif else rc3fits = rc3file
          rc3 = mrdfits(rc3fits,1,h,/silent)
          spawn, ['gzip '+rc3fits], /sh
       end
       else: begin
          print, 'Multiple RC3 FITS files found.'
          return
       end
    endcase

; re-arrange the RC3 catalog to include only useful entries 

    nobj = n_elements(rc3)

; replace all "zeros" (no data) with -999.0.  ignore the first 18
; tags, which correspond to coordinates

    off = 18L

    tags = tag_names(rc3)
    for j = off, n_elements(tags)-off-1L do begin
       if (size(rc3[0].(j),/type) ne 7L) and (tags[j] ne 'PA') and $
         (tags[j] ne 'T') then begin
          w = where(rc3.(j) eq float(0),nw)
          if nw ne 0L then rc3[w].(j) = -999.0
       endif
    endfor

    rc3out = struct_trimtags(rc3,select=['NAME','ALTNAME','PGC','TYPE','TYPESR','T','R25',$
      'PA','BT','BT_CODE','E_BT','BMAG','E_BMAG','M21','E_M21','B_VT','E_B_VT','U_BT','E_U_BT',$
      'HI','W20','E_W20','W50','E_W50','V21'])
    rc3out.name = strcompress(rc3out.name,/remove)

    coords = {RA:  0D, DEC: 0D}
    coords = replicate(coords,nobj)

    info = {D25_MAJ: -999.0, D25_MIN: -999.0, INCLINATION:  -999.0}
    info = replicate(info,nobj)
    
    if not keyword_set(silent) then splog, 'Parsing the RC3 galaxy catalog.'
    for k = 0L, nobj-1L do begin

; parse the coordinates
    
       if rc3[k].ras lt 10.0 then ras = '0'+string(rc3[k].ras,format='(F3.1)') else $
         ras = string(rc3[k].ras,format='(F4.1)')
       rahms = string(rc3[k].rah,format='(I2.2)')+':'+string(rc3[k].ram,format='(I2.2)')+':'+ras
       dechms = rc3[k].de_+string(rc3[k].ded,format='(I2.2)')+':'+$
         string(rc3[k].dem,format='(I2.2)')+':'+string(rc3[k].des,format='(I2.2)')

       coords[k].ra = 15D*hms2dec(rahms)
       coords[k].dec = hms2dec(dechms)

; other quantities
       if (rc3[k].r25 gt -900.0) then begin
       
          ratio = 1.0/10.0^(rc3[k].r25) ; "b" divided by "a"
          quantity = sqrt((1.0-ratio^2.0)/0.96)
          if quantity lt 1.0 then info[k].inclination = asin(quantity<1)*!radeg ; inclination angle [degrees]

       endif

       if (rc3[k].r25 gt -900.0) and (rc3[k].d25 gt -900.0) then begin

          info[k].d25_maj = 0.1 * 10^rc3[k].d25     ; major axis diameter [arcmin]
          info[k].d25_min = ratio * info[k].d25_maj ; minor axis diameter [arcmin]

       endif
       
    endfor

; according to NED the RC3 coordinates are incorrect for UGC05491

    match = where(strmatch(rc3out.altname,'*5491*') eq 1B,nmatch)
    coords[match].ra = 15D*hms2dec('10:11:58.2')
    coords[match].dec = hms2dec('+58:51:52')
    
    rc3out = struct_addtags(coords,rc3out)
    rc3out = struct_addtags(rc3out,info)

    outname = 'rc3_parsed.fits'
    splog, 'Writing '+outname+'.'
    mwrfits, rc3out, outname, /create
    spawn, ['gzip -f '+outname], /sh
    
    popd

return
end

