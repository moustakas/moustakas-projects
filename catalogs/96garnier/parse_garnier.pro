pro parse_garnier, silent=silent
; jm03may4uofa

    garnierpath = filepath('',root_dir=getenv('CATALOGS_DIR'),subdirectory='garnier')
    pushd, garnierpath

    garnierfile = findfile('garnier.fits*',count=fcount)
    if not keyword_set(silent) then splog, 'Reading the original Garnier galaxy catalog.'
    case fcount of
       0L: begin
          print, 'GARNIER FITS file not found.'
          return
       end
       1L: begin
          garnierfile = garnierfile[0]
          if strmatch(garnierfile,'*.gz') eq 1B then begin
             spawn, ['gunzip '+garnierfile], /sh
             garnierfits = strmid(garnierfile,0,strpos(garnierfile,'.gz'))
          endif else garnierfits = garnierfile
          garnier = mrdfits(garnierfits,1,h,/silent)
          spawn, ['gzip '+garnierfits], /sh
       end
       else: begin
          print, 'Multiple GARNIER FITS files found.'
          return
       end
    endcase

; re-arrange the Garnier catalog to include only useful entries.
; precess the coordinates to J2000 

    nobj = n_elements(garnier)

; replace all "zeros" (no data) with -999.0.  ignore the first 9 tags,
; which correspond to coordinates or galaxy names 

    off = 9L
    
    tags = tag_names(garnier)
    for j = off, n_elements(tags)-off-1L do begin
       if (size(garnier[0].(j),/type) ne 7L) and (tags[j] ne 'PA') then begin
          w = where(garnier.(j) eq float(0),nw)
          if nw ne 0L then garnier[w].(j) = -999.0
       endif
    endfor

    garnierout = struct_trimtags(garnier,select=['NAME','PA','BT','E_BT','LOGR25'])

; parse the coordinates and other quantities
    
    coords = {RA:  '', DEC: ''}
    coords = replicate(coords,nobj)

    info = {D25_MAJ: -999.0, D25_MIN: -999.0, INCLINATION:  -999.0}
    info = replicate(info,nobj)
    
    if not keyword_set(silent) then splog, 'Parsing the Garnier galaxy catalog.'
    for k = 0L, nobj-1L do begin

; coordinates
       
       if garnier[k].ras lt 10.0 then ras = '0'+string(garnier[k].ras,format='(F3.1)') else $
         ras = string(garnier[k].ras,format='(F4.1)')
       ra = string(garnier[k].rah,format='(I2.2)')+':'+string(garnier[k].ram,format='(I2.2)')+':'+ras

       dec = garnier[k].de_+string(garnier[k].ded,format='(I2.2)')+':'+$
         string(garnier[k].dem,format='(I2.2)')+':'+string(garnier[k].des,format='(I2.2)')

       ra2000 = 15.0*im_hms2dec(ra)
       dec2000 = im_hms2dec(dec)
       precess, ra2000, dec2000, 1950.0, 2000.0

       coords[k].ra = strjoin(strsplit(im_dec2hms(ra2000/15.0),' ',/extract),':')
       coords[k].dec = strjoin(strsplit(im_dec2hms(dec2000),' ',/extract),':')

; other quantities

       if (garnier[k].logr25 gt -900.0) then begin

          ratio = 1.0/10.0^(garnier[k].logr25) ; "b" divided by "a"
          quantity = sqrt((1.0-ratio^2.0)/0.96)
          if quantity lt 1.0 then info[k].inclination = asin(quantity)*!radeg ; inclination angle [degrees]

       endif

       if (garnier[k].logr25 gt -900.0) and (garnier[k].logd25 gt -900.0) then begin
          
          info[k].d25_maj = 0.1 * 10^garnier[k].logd25 ; major axis diameter [arcmin]
          info[k].d25_min = ratio * info[k].d25_maj    ; minor axis diameter [arcmin]

       endif
          
    endfor

; append everything

    garnierout = struct_addtags(coords,garnierout)
    garnierout = struct_addtags(garnierout,info)

    outname = 'garnier_parsed.fits'
    splog, 'Writing '+outname+'.'
    mwrfits, garnierout, outname, /create
    spawn, ['gzip '+outname], /sh
    
    popd
    
return
end
