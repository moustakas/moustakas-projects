function read_88tully
; jm02may21uofa

    tully_path = filepath('',root_dir=getenv('CATALOGS_DIR'),subdirectory='88tully')

    pushd, tully_path

    tullyfile = file_search('88tully.fits*',count=fcount)
    case fcount of
       0L: begin
          print, 'TULLY catalog not found.'
          return, -1
       end
       1L: begin
          tullyfile = tullyfile[0]
          if strmatch(tullyfile,'*.gz') eq 1B then begin
             spawn, ['gunzip '+tullyfile], /sh
             tullyfits = strmid(tullyfile,0,strpos(tullyfile,'.gz'))
          endif else tullyfits = tullyfile
          tully = mrdfits(tullyfits,1,/silent)
          spawn, ['gzip '+tullyfits], /sh
       end
       else: begin
          print, 'Multiple TULLY FITS files found.'
          return, -1
       end
    endcase

    popd
    
return, tully
end
