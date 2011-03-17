pro highzea_extract, wfits=wfits, noplot=noplot
; jm06jun28uofa - written

    datapath = highzea_path(/spec2d)
    outpath = highzea_path(/spec1d)
    repeatpath = highzea_path(/spec1d)+'repeaters/'

    aperture = 10.0

    t0 = systime(1)

; ###########################################################################
    ispec, 'J0811+4716.fits', /deredden, aperture=aperture, optimal=1, opt_check=0, opt_gauss=0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask
; ###########################################################################
    ispec, 'J0824+5032_1.fits', /deredden, aperture=aperture, optimal=1, opt_check=0, opt_gauss=0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask
    ispec, 'J0824+5032_2.fits', /deredden, aperture=aperture, optimal=1, opt_check=0, opt_gauss=0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask
    ispec, 'J0824+5032_3.fits', /deredden, aperture=aperture, optimal=1, opt_check=0, opt_gauss=0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask
; ###########################################################################
    ispec, 'J0826+4305_1.fits', /deredden, aperture=aperture, optimal=1, opt_check=0, opt_gauss=0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask
    ispec, 'J0826+4305_2.fits', /deredden, aperture=aperture, optimal=1, opt_check=0, opt_gauss=0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask
    ispec, 'J0826+4305_3.fits', /deredden, aperture=aperture, optimal=1, opt_check=0, opt_gauss=0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask
    ispec, 'J0826+4305_4.fits', /deredden, aperture=aperture, optimal=1, opt_check=0, opt_gauss=0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask
; ###########################################################################
    ispec, 'J0827+2954_1.fits', /deredden, aperture=aperture, optimal=1, opt_check=0, opt_gauss=0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask
    ispec, 'J0827+2954_2.fits', /deredden, aperture=aperture, optimal=1, opt_check=0, opt_gauss=0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask
    ispec, 'J0827+2954_3.fits', /deredden, aperture=aperture, optimal=1, opt_check=0, opt_gauss=0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask
    ispec, 'J0827+2954_4.fits', /deredden, aperture=aperture, optimal=1, opt_check=0, opt_gauss=0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask
; ###########################################################################
;   ispec, 'J0828+0336.fits', /deredden, aperture=aperture, optimal=1, opt_check=0, opt_gauss=0, $
;     /meanprofile, tracespec=1, traceorder=1, s1, $
;     datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask
; ###########################################################################
    ispec, 'J0944+0930_1.fits', /deredden, aperture=aperture, optimal=1, opt_check=0, opt_gauss=0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask
    ispec, 'J0944+0930_2.fits', /deredden, aperture=aperture, optimal=1, opt_check=0, opt_gauss=0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask
; ###########################################################################
    ispec, 'J1039+4537_1.fits', /deredden, aperture=aperture, optimal=1, opt_check=0, opt_gauss=0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask
    ispec, 'J1039+4537_2.fits', /deredden, aperture=aperture, optimal=1, opt_check=0, opt_gauss=0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask
    ispec, 'J1039+4537_3.fits', /deredden, aperture=aperture, optimal=1, opt_check=0, opt_gauss=0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask
; ###########################################################################
    ispec, 'J1104+5946.fits', /deredden, aperture=aperture, optimal=1, opt_check=0, opt_gauss=0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask
; ###########################################################################
;   ispec, 'J1109-0039.fits', /deredden, aperture=aperture, optimal=1, opt_check=0, opt_gauss=0, $
;     /meanprofile, tracespec=1, traceorder=1, s1, $
;     datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask
; ###########################################################################
    ispec, 'J1125-0145_1.fits', /deredden, aperture=aperture, optimal=1, opt_check=0, opt_gauss=0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask
    ispec, 'J1125-0145_2.fits', /deredden, aperture=aperture, optimal=1, opt_check=0, opt_gauss=0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask
    ispec, 'J1125-0145_3.fits', /deredden, aperture=aperture, optimal=1, opt_check=0, opt_gauss=0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask
; ###########################################################################
    ispec, 'J1142+6037_1.fits', /deredden, aperture=aperture, optimal=1, opt_check=0, opt_gauss=0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask
    ispec, 'J1142+6037_2.fits', /deredden, aperture=aperture, optimal=1, opt_check=0, opt_gauss=0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask
; ###########################################################################
    ispec, 'J1248+0601.fits', /deredden, aperture=aperture, optimal=1, opt_check=0, opt_gauss=0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask
; ###########################################################################
    ispec, 'J1359+5137.fits', /deredden, aperture=aperture, optimal=1, opt_check=0, opt_gauss=0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask
; ###########################################################################
    ispec, 'J1506+5402_1.fits', /deredden, aperture=aperture, optimal=1, opt_check=0, opt_gauss=0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask
    ispec, 'J1506+5402_2.fits', /deredden, aperture=aperture, optimal=1, opt_check=0, opt_gauss=0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask
    ispec, 'J1506+5402_3.fits', /deredden, aperture=aperture, optimal=1, opt_check=0, opt_gauss=0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask
; ###########################################################################
    ispec, 'J1604+3939_1.fits', /deredden, aperture=aperture, optimal=1, opt_check=0, opt_gauss=0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask
    ispec, 'J1604+3939_2.fits', /deredden, aperture=aperture, optimal=1, opt_check=0, opt_gauss=0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask
; ###########################################################################
    ispec, 'J1634+4619.fits', /deredden, aperture=aperture, optimal=1, opt_check=0, opt_gauss=0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask
; ###########################################################################
    ispec, 'J1635+4709_1.fits', /deredden, aperture=aperture, optimal=1, opt_check=0, opt_gauss=0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask
    ispec, 'J1635+4709_2.fits', /deredden, aperture=aperture, optimal=1, opt_check=0, opt_gauss=0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask
; ###########################################################################
    ispec, 'J1713+2817.fits', /deredden, aperture=aperture, optimal=1, opt_check=0, opt_gauss=0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask
; ###########################################################################
    ispec, 'J2140+1209.fits', /deredden, aperture=aperture, optimal=1, opt_check=0, opt_gauss=0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask
; ###########################################################################

    splog, 'Total time to extract all HIGHZEA spectra (minutes):', (systime(1)-t0)/60.0

return
end
