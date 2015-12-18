pro extract_eva, wfits=wfits, noplot=noplot
; jm04mar18uofa

    outpath = cwd();+'eva/'
    
    ispec, 'fwcsra.0446.fits.gz', /deredden, /heliocor, aperture=2.5, $
      /meanprofile, /tracespec, outpath=outpath, outname='ngc_4441_nuclear_002.5.ms.fits', $
      noplot=noplot, wfits=wfits, /gzip, /telluric
    ispec, 'fwcsra.0447_2.fits.gz', /deredden, /heliocor, aperture=50, $
      /meanprofile, /tracespec, outpath=outpath, outname='ngc_4441_drift_090_050.ms.fits', $
      noplot=noplot, wfits=wfits, /gzip, /telluric

    ispec, 'fwcsra.0450.fits.gz', /deredden, /heliocor, aperture=2.5, $
      /meanprofile, /tracespec, outpath=outpath, outname='ugc_08264_nuclear_002.5.ms.fits', $
      noplot=noplot, wfits=wfits, /gzip, /telluric
    ispec, 'fwcsra.0451_2.fits.gz', /deredden, /heliocor, aperture=50, $
      /meanprofile, /tracespec, outpath=outpath, outname='ugc_08264_drift_040_050.ms.fits', $
      noplot=noplot, wfits=wfits, /gzip, /telluric

    ispec, 'fwcsra.0454.fits.gz', /deredden, /heliocor, aperture=2.5, $
      /meanprofile, /tracespec, outpath=outpath, outname='ngc_4194_nuclear_002.5.ms.fits', $
      noplot=noplot, wfits=wfits, /gzip, /telluric
    ispec, 'fwcsra.0455_2.fits.gz', /deredden, /heliocor, aperture=50, $
      /meanprofile, /tracespec, outpath=outpath, outname='ngc_4194_drift_096_050.ms.fits', $
      noplot=noplot, wfits=wfits, /gzip, /telluric

    ispec, 'fwcsra.0458.fits.gz', /deredden, /heliocor, aperture=2.5, $
      /meanprofile, /tracespec, outpath=outpath, outname='ngc_5607_nuclear_002.5.ms.fits', $
      noplot=noplot, wfits=wfits, /gzip, /telluric
    ispec, 'fwcsra.0459_2.fits.gz', /deredden, /heliocor, aperture=50, $
      /meanprofile, /tracespec, outpath=outpath, outname='ngc_5607_drift_060_050.ms.fits', $
      noplot=noplot, wfits=wfits, /gzip, /telluric

; write postscript and generate the web page

    html_outpath = '/home/ioannis/public_html/eva/'
    plot1dspec, file_search(['*drift*','*nuclear*']), postscript=wfits, outpath=html_outpath

    htmlbase = 'eva'
    html_path = '/home/ioannis/public_html/'
    if keyword_set(wfits) then im_ps2html, htmlbase, html_path=html_path, _extra=extra
    
return
end
