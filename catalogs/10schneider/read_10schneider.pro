function read_10schneider, silent=silent
; jm10jun15ucsd - 
    return, mrdfits(getenv('CATALOGS_DIR')+$
      '/10schneider/10schneider.fits.gz',1,/silent)
end
