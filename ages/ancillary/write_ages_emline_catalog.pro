pro write_ages_emline_catalog
; jm08aug25nyu - generate an emission-line catalog 

    scale = 1D17
    snrcut = 3.0
    chi2cut = 5.0
    sigmacut = 500.0
    
    vv = 'v1.0' ; uses v2.0 emission-line catalog
    
    version = ages_version(/ispec)
    apath = ages_path(/analysis)
    outpath = ages_path(/catalogs)
    
    aa = read_ages(/ancillary)
    ii = read_ages(/ispec)
;   ii = read_ages(/ispec,/notweak)
    uu = read_ages(/unfluxed)
    cc = mrdfits(apath+'catalog.cat.noguidestars.fits.gz',1)
    nobj = n_elements(cc)

    spherematch, aa.ra, aa.dec, 15.0D*cc.cat_ra, $
      cc.cat_dec, 1.5/3600.0, m1, m2
;   niceprint, aa[m1[0:10]].ra, cc[m2[0:10]].cat_ra*15.0D

    lines = ['OII_3727','H_BETA','OIII_5007','H_ALPHA','NII_6584']
    ewlines = lines+'_EW'
    nline = n_elements(lines)

    moretags = replicate({ra: -999.0D, dec: -999.0D, z: -999.0, $
      pass: -999, aper: -999},nobj)
    nmore = n_tags(moretags)

    flux = struct_addtags(struct_addtags(moretags,$
      mrd_struct(lines,replicate('[-999.0,-999.0]',nobj),nobj)),$
      mrd_struct(ewlines,replicate('[-999.0,-999.0]',nobj),nobj))
    format = strarr(2,nmore+2*nline)
    format[*,0:nmore-1L] = [['RA','F12.7'],['DEC','F12.7'],$ ; for the output file
      ['Z','F12.5'],['PASS','I5'],['APER','I5']]
    for iline = 0L, nline-1L do format[*,nmore+iline] = $ ; flux
      [transpose([[lines[iline]],['2F12.5']])] 
    for iline = 0L, nline-1L do format[*,nmore+nline+iline] = $ ; EWs
      [transpose([[ewlines[iline]],['2F12.5']])]

    rawflux = ii
    
    alltags = tag_names(ii)
    for iline = 0L, nline-1L do begin

       fluxtag = where(lines[iline] eq alltags)
       ewtag = where(lines[iline]+'_EW' eq alltags)
       chi2tag = where(lines[iline]+'_CHI2' eq alltags)
       sigmatag = where(lines[iline]+'_SIGMA' eq alltags)

       fluxcrap = where($
         (aa.pass eq 106) or (aa.pass eq 110) or $
         (aa.pass eq 209) or (aa.pass eq 310) or $
         (aa.pass eq 311) or $
         ((ii.(fluxtag))[0,*]/(ii.(fluxtag))[1,*] lt snrcut) or $
         (ii.(chi2tag) gt chi2cut) or (ii.(sigmatag) ge sigmacut),ncrap,comp=fluxgood)

       rawflux[fluxcrap].(fluxtag)[0] = -999.0
       rawflux[fluxcrap].(fluxtag)[1] = -999.0

       rawflux[fluxcrap].(ewtag)[0] = -999.0
       rawflux[fluxcrap].(ewtag)[1] = -999.0

; scale the fluxes       
       
       rawflux[fluxgood].(fluxtag)[0] = scale*rawflux[fluxgood].(fluxtag)[0]
       rawflux[fluxgood].(fluxtag)[1] = scale*rawflux[fluxgood].(fluxtag)[1]

    endfor

; copy the parsed results to the large output structure
    
    junk = flux[m2] & struct_assign, rawflux[m1], junk, /nozero & flux[m2] = junk
    flux.ra = cc.cat_ra & flux.dec = cc.cat_dec

;   niceprint, ii[m1[0:10]].ra, flux[m2[0:10]].ra*15.0D, ii[m1[0:10]].z, flux[m2[0:10]].z

; write out

    outfile = outpath+'catalog.moustakas.emline.'+vv+'.fits'
    splog, 'Writing '+outfile
    mwrfits, flux, outfile, /create
    
    outfile = outpath+'catalog.moustakas.emline.'+vv
    splog, 'Writing '+outfile

    openw, lun, outfile, /get_lun
    printf, lun, '# J. Moustakas, NYU'
    printf, lun, '# Contact john.moustakas@gmail.com with questions'
    printf, lun, '# '+im_today()
    printf, lun, '# '
    printf, lun, '# Catalog of strong emission-line fluxes and equivalent widths based '
    printf, lun, '# on the '+version+' AGES emission-line catalog.  The first five columns '
    printf, lun, '# identify the AGES spectrum used and the adopted redshift. The next '
    printf, lun, '# '+string(2*nline,format='(I0)')+' columns give the fluxes and errors (in units of 1E-17 erg/s/cm2) '
    printf, lun, '# of the strong emission lines based on simultaneous multi-Gaussian fits '
    printf, lun, '# to the continuum-subtracted spectra.  Finally, the last '+string(2*nline,format='(I0)')+' columns '
    printf, lun, '# give the corresponding rest-frame equivalent width (EW) of each '
    printf, lun, '# line in Angstroms.  Null values are indicated with -999.'
    printf, lun, '# '
    printf, lun, '# A line had to satisfy the following criteria to be included:'
    printf, lun, '#   * S/N>3'
    printf, lun, '#   * Reduced chi^2<5'
    printf, lun, '#   * Gaussian sigma width <500 km/s'
    printf, lun, '# '
    printf, lun, '# Furthermore, measurements from the following plates were excluded because '
    printf, lun, '# they were not spectrophotometrically calibrated: 106, 110, 209, 310, and '
    printf, lun, '# 311.  Several other plates have spectrophotometric calibration issues, but '
    printf, lun, '# have been included in the spirit of being complete; contact me for details.  '
    printf, lun, '# I have also measured EWs from the *unfluxed* spectra; again, contact me if '
    printf, lun, '# you need these.'
    printf, lun, '# '
    printf, lun, '# Note that these quality cuts do not remove *all* the spurious line-'
    printf, lun, '# measurements; **visual inspection is strongly recommended!**  In '
    printf, lun, '# particular, spectra with fewer than 2-3 well-detected lines may not '
    printf, lun, '# be bona fide emission-line galaxies.'
    printf, lun, '# '
    printf, lun, '# Future versions of this catalog may include additional lines and possibly '
    printf, lun, '# upper limits on the undetected lines.'
    printf, lun, '# '
    printf, lun, '#  1 RA'
    printf, lun, '#  2 DEC'
    printf, lun, '#  3 Z'
    printf, lun, '#  4 PASS'
    printf, lun, '#  5 APER'
    counter = 6
    for iline = 0L, nline-1L do begin
       printf, lun, '# '+string(counter,format='(I2)')+' '+lines[iline]+'_FLUX'
       counter = counter+1
       printf, lun, '# '+string(counter,format='(I2)')+' '+lines[iline]+'_FLUX_ERR'
       counter = counter+1
    endfor
    for iline = 0L, nline-1L do begin
       printf, lun, '# '+string(counter,format='(I2)')+' '+lines[iline]+'_EW'
       counter = counter+1
       printf, lun, '# '+string(counter,format='(I2)')+' '+lines[iline]+'_EW_ERR'
       counter = counter+1
    endfor
    struct_print, flux, format=format, lun=lun, /no_head
    free_lun, lun
;   for jj = 0L, nobj-1L do printf, lun, 
;   README

stop    
    

return
end
    
