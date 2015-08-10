pro decals_cosmos_cfhtls
; jm15aug08siena - match the DECaLS photometry of the COSMOS field to the
; CFHTLS/D2 photometry from Gwyn+

    path = getenv('DESI_ROOT')+'/target/analysis/deep2/v2.0/'
    cfhtlsfile = path+'cfhtls-D2.fits'

    splog, 'Reading '+cfhtlsfile
    cat = mrdfits(cfhtlsfile,1)
    ngal = n_elements(cat)

    dr1mdir = getenv('DECALS_DIR_DR1')+'m/'
    allbricks = mrdfits(dr1mdir+'dr1m-bricks.fits',1)

    spherematch, cat.ra, cat.dec, allbricks.ra, allbricks.dec, 0.25D, m1, m2, maxmatch=0

    keep = m2[uniq(m2,sort(m2))]
    bricks = allbricks[keep]

    catfiles = dr1mdir+'tractor/'+strmid(bricks.brickname,0,3)+'/tractor-'+$
      bricks.brickname+'.fits'
       
    for ib = 0, n_elements(bricks)-1 do begin
       t1 = mrdfits(catfiles[ib],1)
       if n_elements(tractor) eq 0L then tractor = t1 else $
         tractor = [t1,tractor]
    endfor
    
    spherematch, tractor.ra, tractor.dec, cat.ra, cat.dec, 1D/3600.0, m1, m2
    out = im_empty_structure(tractor[0],ncopies=ngal,empty_value=-999)
    out[m2] = im_struct_assign(tractor[m1],out[m2])

    splog, 'Writing '+path+'decals-cfhtls-D2.fits'
    mwrfits, out, path+'decals-cfhtls-D2.fits', /create
    
return
end
