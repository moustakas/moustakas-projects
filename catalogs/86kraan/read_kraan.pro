function read_kraan
; jm02dec10uofa

    kraanpath = filepath('',root_dir=getenv('CATALOGS_DIR'),subdirectory='kraan')

    pushd, kraanpath

    kraanfile = findfile('kraan.fits*',count=fcount)
    case fcount of
       0L: begin
          print, 'KRAAN-KORTEWEG catalog not found.'
          return, -1
       end
       1L: begin
          kraanfile = kraanfile[0]
          if strmatch(kraanfile,'*.gz') eq 1B then begin
             spawn, ['gunzip '+kraanfile], /sh
             kraanfits = strmid(kraanfile,0,strpos(kraanfile,'.gz'))
          endif else kraanfits = kraanfile
          kraan = mrdfits(kraanfits,1,/silent)
          spawn, ['gzip '+kraanfits], /sh
       end
       else: begin
          print, 'Multiple KRAAN-KORTEWEG FITS files found.'
          return, -1
       end
    endcase

    popd
    
return, kraan
end
