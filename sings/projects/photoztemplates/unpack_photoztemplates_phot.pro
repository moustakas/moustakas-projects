pro unpack_photoztemplates_phot

    path = sings_path(/proj)+'photoztemplates/'
    ff = file_search('NGC*.phot',count=nff)

    phot = {galaxy: '', z: 0.0, maggies: fltarr(17), ivarmaggies: fltarr(17)}
    phot = replicate(phot,nff)
    phot.galaxy = repstr(ff,'.phot','')
    
    for ii = 0, nff-1 do begin
       all = djs_readilines(ff[ii])
       nrow = n_elements(all)
       mag = fltarr(nrow)
       for jj = 0, nrow-1 do mag[jj] = float((strsplit(all[jj],' ',/extract))[1])

       phot[ii].maggies[0:nrow-1] = mag2maggies(mag,magerr=mag*0+0.05,ivarmaggies=ivar)
       phot[ii].ivarmaggies[0:nrow-1] = ivar
    endfor

; get the redshifts and also write out the optical spectra
;   info = read_sings_gandalf(/drift56)
    n4214 = where(phot.galaxy eq 'NGC4214')
    phot[n4214].z = 0.000970000 ; special case
    
    d56 = mrdfits(path+'sings_drift56_models.fits.gz',1)
    match, strtrim(d56.galaxy,2), strtrim(phot.galaxy,2), m1, m2
    phot[m2].z = d56.zabs[m1]

; normalize the optical spectrum to the g-band flux and convert to mAB 
    filters = photoztemplates_filterlist()
    normfilt = 'sdss_g0.par'
    this = (where(normfilt eq filters))[0]

    wave = d56.wave
    flux = fltarr(n_elements(wave),nff)
    isdata = intarr(n_elements(wave),nff)

    flux[*,m2] = d56.flux[*,m1]
    isdata[*,m2] = d56.isdata[*,m1]

    for ii = 0, n_elements(m2)-1 do begin
       optmaggies = (k_project_filters(k_lambda_to_edges(wave),flux[*,m2[ii]],$
         filterlist=normfilt))[0]
       norm = phot[m2[ii]].maggies[this]/optmaggies
       splog, phot[m2[ii]].galaxy, norm
       flux[*,m2[ii]] = flux[*,m2[ii]]*norm
; convert to mab
       flux[*,m2[ii]] = flux[*,m2[ii]]*wave^2/im_light(/ang)
       flux[*,m2[ii]] = -2.5*alog10(flux[*,m2[ii]]>1D-50)-48.6 ; [AB mag]
    endfor
    spec = {galaxy: phot.galaxy, zabs: phot.z, wave: wave, $
      flux: flux, isdata: isdata}

    im_mwrfits, phot, 'phot.fits', /clobber
    im_mwrfits, spec, 'optspec.fits', /clobber
    
return
end
    
