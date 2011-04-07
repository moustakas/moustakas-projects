function read_bcgs_photometry, sdss=sdss, ivarmaggies=ivarmaggies, $
  zobj=zobj, galaxy=galaxy, filterlist=filterlist
; jm11apr06ucsd - simple wrapper to pull out the photometry we need

    bcgspath = ages_path(/proj)+'bcgs/'
    
    if keyword_set(sdss) then begin
       phot = mrdfits(bcgspath+'bcgs_sdss_phot.fits.gz',1)
       maggies = phot.maggies
       ivarmaggies = phot.ivarmaggies
       zobj = phot.z
       galaxy = phot.galaxy
       filterlist = [sdss_filterlist(),twomass_filterlist()]
    endif else begin
; do not use ch4 when it shifts redward of 3 microns, and do not use
; the Bw-band blueward of ~2200 A
       vv = 'v3'
       phot = mrdfits(bcgspath+'bcgs_photometry_'+vv+'.fits.gz',1)
       sample = rsex(bcgspath+'bcgs_sample_'+vv+'.sex')
       zobj = sample.z
       bootes_to_maggies, phot, maggies, ivarmaggies, $
         use_aper='08', filterlist=filterlist
       galaxy = hogg_iau_name(sample.ra,sample.dec,'Bootes')

       keep = where((strmatch(filterlist,'*ufilter*') eq 0) and $
         (strmatch(filterlist,'*bok*') eq 0))
       filterlist = filterlist[keep]
       maggies = maggies[keep,*]
       ivarmaggies = ivarmaggies[keep,*]
       
       Bw = (where(strmatch(filterlist,'*Bw*',/fold)))[0]
       Bwtoss = where(zobj gt 1.0)

       ch3 = (where(strmatch(filterlist,'*ch3*')))[0]
       ch4 = (where(strmatch(filterlist,'*ch4*')))[0]
       ivarmaggies[Bw,Bwtoss] = 0.0
       ivarmaggies[ch3,*] = 0.0
       ivarmaggies[ch4,*] = 0.0
    endelse

return, maggies
end
