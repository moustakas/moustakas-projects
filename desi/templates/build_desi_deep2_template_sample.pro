function get_ugriz, cat
; compute dust-corrected ugriz magnitudes
    deep2_to_maggies, cat, mm, ii
    offset = 3 ; offset from the BRI photometry
    ugriz = fltarr(5,n_elements(cat))-99.0
    for ii = 0, 4 do begin
       good = where(mm[offset+ii,*] gt 0)
       ugriz[ii,good] = -2.5*alog10(mm[offset+ii,good])
    endfor
return, ugriz
end

pro build_desi_deep2_template_sample
; jm13dec18siena - build the sample of DEEP2 galaxies we will use to
;   construct DESI templates
; jm14mar11siena - major update

    templatepath = getenv('IM_PROJECTS_DIR')+'/desi/templates/'

    snrcut = 3.0
    
    zcat = read_deep2_zcat(photo=photo)
    line = read_deep2(/ppxf)
    linefixoii = read_deep2(/ppxf,/fixoii)
    kised = mrdfits(templatepath+'desi_deep2_fsps_v2.4_miles_'+$
      'chab_charlot_sfhgrid01_kcorr.z0.0.fits.gz',1)

    ugriz = get_ugriz(photo)

    
    
    these = where($
      zcat.zbest gt 0.7 and zcat.zbest lt 1.6 and $
      total(maggies gt 0,1) ge 3 and $
      -2.5*alog10(maggies[1,*]) lt 24.0,ngal)
;     ispec.oii_3727_1_ew[0] gt 1.0 and $
;     ispec.oii_3727_2_ew[0] gt 1.0 and $
;     ispec.oii_3727_1[0]/ispec.oii_3727_1[1] gt 5.0 and $
;     ispec.oii_3727_2[0]/ispec.oii_3727_2[1] gt 5.0,ngal)
    splog, 'Sample ', ngal
    
    
stop    
    
    
    kcorr = read_deep2(/kcorr)
    ngal = n_elements(zcat)
    
    
    deep2_to_maggies, phot, maggies, ivarmaggies, filterlist=filt
    nfilt = n_elements(filt)

    zcat = struct_addtags(zcat,replicate({maggies: fltarr(nfilt), $
      ivarmaggies: fltarr(nfilt), source: ''},n_elements(zcat)))
    zcat.maggies = maggies
    zcat.ivarmaggies = ivarmaggies
    zcat.source = phot.source


    these = where($
      zcat.zbest gt 0.7 and zcat.zbest lt 1.6 and $
      total(maggies gt 0,1) ge 3 and $
      -2.5*alog10(maggies[1,*]) lt 24.0,ngal)
;     ispec.oii_3727_1_ew[0] gt 1.0 and $
;     ispec.oii_3727_2_ew[0] gt 1.0 and $
;     ispec.oii_3727_1[0]/ispec.oii_3727_1[1] gt 5.0 and $
;     ispec.oii_3727_2[0]/ispec.oii_3727_2[1] gt 5.0,ngal)
    splog, 'Sample ', ngal

; build a poor-mans desi sample by applying on [OII] flux simple cut
; >5e-17 cgs 
;* convert the fluxes to physical units
;* build the full header - add an EXTNAME tag
;* document each tag of the FITS table
;* more tags - integrated [OII] flux, flux in each component, sigma, EW, continuum flux

    im_mwrfits, zcat[these], templatepath+'deep2_zcat.fits', /clobber
    im_mwrfits, kcorr[these], templatepath+'deep2_kcorr.fits', /clobber
    im_mwrfits, ispec[these], templatepath+'deep2_ispec.fits', /clobber
    
return
end
