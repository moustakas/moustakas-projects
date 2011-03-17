pro atlas_literature_photometry, write=write
; jm05mar02uofa
; compile B-band photometry from the literature

    if keyword_set(write) then postscript = 1L
    
    light = 2.99792458D5 ; speed of light [km/s]
    
    nedpath = atlas_path(/ned)
    outpath = atlas_path(/analysis)

    atlas = mrdfits(outpath+'atlas_ancillary_data.fits.gz',1,/silent) ; this must match LEDA!
    ngalaxy = n_elements(atlas)

    ra = 15.0*im_hms2dec(atlas.ra)
    dec = im_hms2dec(atlas.dec)

; ---------------------------------------------------------------------------
; initialize the output data structure
; ---------------------------------------------------------------------------

    photdata = {$
      galaxy:       '', $
      ra:           '', $
      dec:          '', $
      rc3_B:     -999.0,$
      rc3_B_err: -999.0,$
      B:         -999.0,$
      B_err:     -999.0,$
      B_ref:        ''}
    photdata = replicate(photdata,ngalaxy)

    photdata.galaxy = strtrim(atlas.galaxy,2)
    photdata.ra = atlas.ra
    photdata.dec = atlas.dec
    photdata.rc3_b = atlas.rc3_b
    photdata.rc3_b_err = atlas.rc3_b_err

; ---------------------------------------------------------------------------
; read the Tully et al. (1996) - Ursa Major - catalog
; ---------------------------------------------------------------------------

    nophot = where(photdata.B lt -900.0,nnophot)
    
    print, format='("Reading the Ursa Major data . . . ",$)'
    umajor = read_96tully()

    ura = 15.0*im_hms2dec(umajor.ra)
    ude = im_hms2dec(umajor.dec)

    ntot = im_djs_angle_match(ra[nophot],dec[nophot],ura,ude,dtheta=30.0/3600.0,$
      units='degrees',mcount=mcount,mindx=mindx,mdist=mdist)
    good = where(mindx ne -1L,ngood)
    print, 'found '+string(ngood,format='(I0)')+'/'+string(nnophot,format='(I0)')+' galaxies.'

    srt = sort(mdist[good])
    niceprint, photdata[nophot[good[srt]]].galaxy, strcompress(umajor[mindx[good[srt]]].galaxy,/remove), $
      mdist[good[srt]]*3600.0
    
    photdata[nophot[good]].B = umajor[mindx[good]].B
    photdata[nophot[good]].B_err = umajor[mindx[good]].B_err
    photdata[nophot[good]].B_ref = 'Tully et al. (1996)'
    print

; ---------------------------------------------------------------------------
; read the Gil de Paz (2003) catalog
; ---------------------------------------------------------------------------

    nophot = where(photdata.B lt -900.0,nnophot)

    print, format='("Reading the Gil de Paz data . . . ",$)'
    gil = read_03gildepaz()

    ura = 15.0*im_hms2dec(gil.ra)
    ude = im_hms2dec(gil.dec)

    ntot = im_djs_angle_match(ra[nophot],dec[nophot],ura,ude,dtheta=30.0/3600.0,$
      units='degrees',mcount=mcount,mindx=mindx,mdist=mdist)
    good = where(mindx ne -1L,ngood)
    print, 'found '+string(ngood,format='(I0)')+'/'+string(nnophot,format='(I0)')+' galaxies.'

    srt = sort(mdist[good])
    niceprint, photdata[nophot[good[srt]]].galaxy, strcompress(gil[mindx[good[srt]]].galaxy,/remove), $
      mdist[good[srt]]*3600.0
    print
    
    photdata[nophot[good]].B = gil[mindx[good]].Bmag
    photdata[nophot[good]].B_err = gil[mindx[good]].E_Bmag
    photdata[nophot[good]].B_ref = 'Gil de Paz et al. (2003)'

; ---------------------------------------------------------------------------
; read the de Jong (1994) catalog - don't use this!  poor coordinates!
; ---------------------------------------------------------------------------

;   nophot = where(photdata.B lt -900.0,nnophot)
;
;   print, format='("Reading the de Jong data . . . ",$)'
;   dejong = read_94deJong()
;
;   ura = 15.0*im_hms2dec(dejong.ra)
;   ude = im_hms2dec(dejong.dec)
;
;   ntot = im_djs_angle_match(ra,dec,ura,ude,dtheta=30.0/3600.0,$
;     units='degrees',mcount=mcount,mindx=mindx,mdist=mdist)
;   good = where(mindx ne -1L,ngood)
;   print, 'found '+string(ngood,format='(I0)')+'/'+string(nnophot,format='(I0)')+' galaxies.'
;
;   srt = sort(mdist[good])
;   niceprint, photdata[nophot[good[srt]]].galaxy, dejong[mindx[good[srt]]].galaxy, $
;     dejong[mindx[good[srt]]].alt_galaxy, mdist[good[srt]]*3600.0
;   print
;   
;   photdata[nophot[good]].B = dejong[mindx[good]].MB
;   photdata[nophot[good]].B_err = dejong[mindx[good]].E_MB
;   photdata[nophot[good]].B_ref = 'de Jong et al. (1994)'

; ---------------------------------------------------------------------------
; write out    
; ---------------------------------------------------------------------------
    
    if keyword_set(write) then begin

       outfile = 'atlas_photometry.fits'
       splog, 'Writing '+outpath+outfile
       mwrfits, photdata, outpath+outfile, /create
       spawn, ['gzip -f '+outpath+outfile], /sh

    endif

stop    
    
    good = where(photdata.b gt -900 and photdata.rc3_b gt -900,ngood)
    plot, photdata[good].b, photdata[good].rc3_b, ps=4, xsty=3, ysty=3, $
      xrange=[10,17], yr=[10,17]
    oplot, !x.crange, !y.crange, line=0, thick=2



stop    

return
end    
