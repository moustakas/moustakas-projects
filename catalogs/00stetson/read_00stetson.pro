function read_00stetson
; jm08jul18nyu - read the Stetson (2000) photometric standard star
;                catalog 
    return, mrdfits(getenv('CATALOGS_DIR')+$
      '/00stetson/phstd_stetson.tfits',1,/silent)
end
