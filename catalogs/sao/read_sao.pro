function read_sao, silent=silent
; jm10jun15ucsd - 
    return, mrdfits(getenv('CATALOGS_DIR')+$
      '/sao/sao.fits.gz',1,/silent)
end
