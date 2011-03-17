function read_garnier, silent=silent
; jm03apr29upfa

    garnierpath = filepath('',root_dir=getenv('CATALOGS_DIR'),subdirectory='garnier')
    pushd, garnierpath

    garnierfile = findfile('garnier_parsed.fits*',count=fcount)
    if not keyword_set(silent) then splog, 'Reading the Garnier galaxy catalog.'
    case fcount of
       0L: begin
          print, 'GARNIER FITS file not found.'
          return, -1
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
          return, -1
       end
    endcase

    popd

return, garnier
end
