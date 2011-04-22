function read_hizea_photometry, ivarmaggies=ivarmaggies, zobj=zobj, $
  galaxy=galaxy, filterlist=filterlist
; jm10dec20ucsd - simple wrapper to pull out the photometry we want
; for the HIZEA project; called HIZEA_ISEDFIT 

    filterlist = hizea_filterlist()
    
    hizeapath = getenv('HIZEA_DATA')+'/'
    info = mrdfits(hizeapath+'sdss/hizea_sdss_simplefit_dr7.fits.gz',1)
;   info = mrdfits(hizeapath+'sdss/hizea_spitzer_photofit_v1.0.fit',1)
    sdss = mrdfits(hizeapath+'sdss/hizea_sdss_photo_dr7.fits.gz',1)
    galex = mrdfits(hizeapath+'sdss/hizea_galex_gr6.fits.gz',1)
    ngal = n_elements(info)

    zobj = info.z
    galaxy = hogg_iau_name(sdss.ra,sdss.dec,'SDSS')
    
; line-match everything
    irac1 = mrdfits(hizeapath+'spitzer/hizea_spitzer.fits.gz',1)
    irac = im_empty_structure(irac1,empty_value=-999.0,ncopies=ngal)
    spherematch, info.ra, info.dec, irac1.ra, irac1.dec, 1.0/3600.0, m1, m2
    irac[m1] = irac1[m2]

; convert to maggies; note: the IRAC fluxes are in microJy, so to
; convert to maggies multiply by 10^(-0.4*23.9)
    im_galex_to_maggies, galex, gmaggies, givarmaggies
    sdss_to_maggies, smaggies, sivarmaggies, tsobj=sdss, flux='model'

    imaggies = fltarr(2,ngal)
    iivarmaggies = fltarr(2,ngal)

    ch1good = where(irac.ch1flux gt 0.0)
    ch2good = where(irac.ch2flux gt 0.0)
    imaggies[0,ch1good] = irac[ch1good].ch1flux*10.0^(-0.4*23.9)
    iivarmaggies[0,ch1good] = 1.0/(imaggies[0,ch1good]*0.1)^2
    
    imaggies[1,ch2good] = irac[ch2good].ch2flux*10.0^(-0.4*23.9)
    iivarmaggies[1,ch2good] = 1.0/(imaggies[1,ch2good]*0.1)^2
    
    maggies = [gmaggies,smaggies,imaggies]
    ivarmaggies = [givarmaggies,sivarmaggies,iivarmaggies]

return, maggies
end
