function read_04kara
; jm10mar23ucsd - read the Karachentsev+04 local volume (LV) catalog
    return, mrdfits(getenv('CATALOGS_DIR')+$
      '/04kara/04kara.fits.gz',1,/silent)
end
