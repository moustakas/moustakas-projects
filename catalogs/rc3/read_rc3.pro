function read_rc3, silent=silent
; jm02mar23uofa

    rc3_path = filepath('',root_dir=getenv('CATALOGS_DIR'),subdirectory='rc3')

    rc3file = file_search(rc3_path+'rc3_parsed.fits*',count=fcount)
    if not keyword_set(silent) then splog, 'Reading the RC3 galaxy catalog.'
    case fcount of
       0L: begin
          print, 'RC3 FITS file not found.'
          return, -1
       end
       1L: rc3 = mrdfits(rc3file,1,hh)
       else: begin
          print, 'Multiple RC3 FITS files found.'
          return, -1
       end
    endcase

return, rc3
end
