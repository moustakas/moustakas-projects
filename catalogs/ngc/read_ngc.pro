function read_ngc, silent=silent
; jm15mar19siena - read the NGC catalog

    ngc_path = filepath('',root_dir=getenv('CATALOGS_DIR'),subdirectory='ngc')

    ngcfile = file_search(ngc_path+'ngc.fits*',count=fcount)
    if not keyword_set(silent) then splog, 'Reading the NGC galaxy catalog.'
    case fcount of
       0L: begin
          print, 'NGC FITS file not found.'
          return, -1
       end
       1L: ngc = mrdfits(ngcfile,1,hh)
       else: begin
          print, 'Multiple NGC FITS files found.'
          return, -1
       end
    endcase

; add the RA, DEC coordinates
    ngc = struct_addtags(replicate({ra: 0D, dec: 0D},n_elements(ngc)),ngc)
    ra = ngc.ra1975 & dec = ngc.de1975
    precess, ra, dec, 1975.0, 2000.0
    ngc.ra = ra & ngc.dec = dec

return, ngc
end
