pro stitch_ngc_3344, wfits=wfits
; jm03dec28uofa
; stitch together the East-West scans for NGC 5194 (M51)

; read the spectra; 90 degree slit position angle; East is at lower
; row number

    datapath = cwd()
    fspec1 = 'fscra.2108.fits' ; West
    fspec2 = 'fscra.2109.fits' ; East
    info = iforage([fspec1,fspec2])
    
    spec1 = rd2dspec(fspec1,datapath=datapath,wset=wset)
    spec2 = rd2dspec(fspec2,datapath=datapath)
    nimage = 2L
    
; determine the pixel shift graphically here; shift the western scan
; to match the eastern scan

    ncols = spec1.naxis1
    nrows = spec1.naxis2
    rowaxis = lindgen(nrows)
    bignrows = 2*nrows
    bigindx = lindgen(bignrows)
    
    bigrefspec = bigindx*0.0
    bigspec = bigindx*0.0

    refspec = im_normalize(total(spec2.image,1),rowaxis,/max,normwave=nrows-15L,binsize=14)
    bigrefspec[0:nrows-1L] = refspec

    spec = im_normalize(total(spec1.image,1),rowaxis,/max,normwave=15L,binsize=14)
    bigspec[nrows:bignrows-1L] = spec

    pixshift = -19L ; <-- here

    plot, bigindx, bigrefspec, ps=10, xsty=3, ysty=3, yrange=[0,1.1], xrange=[90,150]
    djs_oplot, bigindx, bigspec, ps=10, color='red'
    djs_oplot, bigindx, shift(bigspec,pixshift), ps=10, color='yellow'
stop
; form the output arrays and co-add the 2D data

    outnrows = bignrows + pixshift

; mean image and sky spectrum; only divide by NIMAGE in the
; overlapping rows
    
    outimage1 = dblarr(ncols,outnrows) & outimage2 = outimage1*0.0
    outimage1[*,0L:nrows-1L] = spec2.image
    outimage2[*,nrows+pixshift:outnrows-1L] = spec1.image

    outskyimage1 = dblarr(ncols,outnrows) & outskyimage2 = outskyimage1*0.0
    outskyimage1[*,0L:nrows-1L] = spec2.sky
    outskyimage2[*,nrows+pixshift:outnrows-1L] = spec1.sky

    denom = outimage1*0.0+1.0
    denom[where((outimage1 gt 0.0) and (outimage2 gt 0.0))] = 2.0
    
    outimage = float(total([ [ [outimage1] ], [ [outimage2] ] ],3)/denom)
    outskyimage = float(total([ [ [outskyimage1] ], [ [outskyimage2] ] ],3)/denom)

    if not keyword_set(wfits) then begin
       
       plot, outimage[600,*], ps=10, xsty=3, ysty=3, yrange=[1E-17,1E-14], thick=2.0
       djs_oplot, outimage1[600,*], ps=10, color='red'
       djs_oplot, outimage2[600,*], ps=10, color='yellow'

    endif

; error map
    
    outsigmap1 = dblarr(ncols,outnrows) & outsigmap2 = outsigmap1*0.0
    outsigmap1[*,0L:nrows-1L] = spec2.sigmamap
    outsigmap2[*,nrows+pixshift:outnrows-1L] = spec1.sigmamap
    outsigmap = float(sqrt(total([ [ [outsigmap1^2.0] ], [ [outsigmap2^2.0] ] ],3)/denom))

; bad pixel mask
    
    outmask1 = fltarr(ncols,outnrows) & outmask2 = outmask1*0.0
    outmask1[*,0L:nrows-1L] = spec2.mask
    outmask2[*,nrows+pixshift:outnrows-1L] = spec1.mask
    outmask = fix(total([ [ [outmask1] ], [ [outmask2] ] ],3))

; telluric spectrum    

    telluric_spec = spec1.telluric_spec
    telluric_head = spec1.telluric_head

; update the header; assume that EXPTIME, DATE-OBS, SCANLEN, PA, CCD
; noise parameters, wavelength parameters, RMAXERR, and RMEANERR, are
; the same; retain the time parameters of the zeroth exposure (UT,
; UTMIDDLE, JD, and ST) and all the ISPEC2D keywords (IRAW, etc.);
; compute the mean AIRMASS, ZD, PARANGLE, RA, and DEC; compute the new
; reference pixel CRPIX2; add a history note
    
    outheader = *spec1.header

    meanra = repstr(strtrim(im_dec2hms(djs_mean(im_hms2dec(info.ra))),2),' ',':')
    meande = repstr(strtrim(im_dec2hms(djs_mean(im_hms2dec(info.dec))),2),' ',':')
    
    sxaddpar, outheader, 'OBJECT', 'N3344 drift'
    sxaddpar, outheader, 'RA', meanra, ' right ascension [HMS]'
    sxaddpar, outheader, 'DEC', meande, ' declination [DMS]'
    sxaddpar, outheader, 'AIRMASS', float(djs_mean(info.airmass)), ' airmass at mid-exposure', format='(F12.4)'
    sxaddpar, outheader, 'ZD', float(djs_mean(info.zd)), after='AIRMASS', $
      ' zenith distance at mid-exposure [degrees]', format='(F12.2)'
    sxaddpar, outheader, 'PARANGLE', float(djs_mean(info.parangle)), after='AIRMASS', $
      ' parallactic angle at mid-exposure [degrees]', format='(F12.2)'
    sxaddpar, outheader, 'CRPIX2', outnrows/2L+1L, ' reference pixel number'
    sxaddhist, "'Images "+strjoin(repstr(info.file,'.fits.gz',''),' and ')+" stitched together "+im_today()+"'", outheader
    
; write out

    if keyword_set(wfits) then begin

       outname = 'fscra.2108-09.fits'
       splog, 'Writing '+datapath+outname+'.'
       wrt2dspec, outname, outimage, outsigmap, outmask, wset=wset, $
         outheader, skyimage=outskyimage, telluric_spec=telluric_spec, $
         telluric_head=telluric_head, datapath=datapath

    endif

    icleanup, spec1
    icleanup, spec2

return
end    
