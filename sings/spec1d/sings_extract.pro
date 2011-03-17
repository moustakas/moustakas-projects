pro sings_extract, wfits=wfits, noplot=noplot
; jm04febuofa
; jm05apr07uofa - the radial strip lengths have been hard-wired for
;                 those objects that have been observed; see
;                 IRS_LL_PARAMS.TXT     
; jm05jul25uofa - updated
; jm06jan26uofa - spectra are now stitched within each observing run
;                 directory; the 55" drift scans extraction diameters
;                 are set to 0.55*R25; holmberg_ii_drift55 (01nov)
;                 removed from further analysis -- totally crappy!
; jm08jan15nyu - additional crappy spectra removed from further
;                analysis; also remove the radial-strip spectrum of
;                NGC3034 from further analysis as it doesn't
;                contain much useful information
; jm08feb06nyu - remove NGC4236/DRIFT56
; jm08oct19nyu - remove NGC2403/NUCLEAR
    
; NGC1482 East/West "knot" have not been extracted    
    
    datapath = sings_path(/spec2d)
    outpath = sings_path(/spec1d)
    repeatpath = sings_path(/spec1d)+'repeaters/'

    nucap = 2.5     ; nuclear aperture [arcsec]
    drift20 = 20.0  ; 20" scan aperture [arcsec]
    drift56 = 56.0  ; 55" scan aperture [arcsec]

    t0 = systime(1)

; ###########################################################################
; no nucleus
; asymmetric spatial profile, centered on the HII Region
    ispec, 'ddo_053_drift_020.fits', /deredden, aperture=drift20, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; D25=92.9; 0.55*D25=51.1
    ispec, 'ddo_053_drift_055.fits', /deredden, aperture=51.1, $
      /meanprofile, tracespec=1, traceorder=1, s2, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; ###########################################################################
; DDO154 - no data; D25=181; 0.55*D25=99.7
; ###########################################################################
; jm08jan15nyu - crappy!
; no nucleus
; D25=208; 0.55*D25=114; crummy
;   ispec, 'ddo_165_drift_056.fits', /deredden, refrow=60, aperture=114.0, $
;     /meanprofile, tracespec=1, traceorder=1, $
;     datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; crummy    
;   ispec, 'ddo_165_drift_020.fits', /deredden, refrow=60, aperture=drift20, $
;     /meanprofile, tracespec=1, traceorder=1, fluxrange=[0,12]*1E-17, s1, $
;     datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; ###########################################################################
; D25=218; 0.55*D25=120; crummy; no nucleus
; jm08jan15nyu - crappy!
;   ispec, 'holmberg_i_drift_020.fits', /deredden, aperture=drift20, $
;     /meanprofile, tracespec=1, traceorder=1, s1, $
;     datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; ###########################################################################
; Holmberg II; no data; D25=477; 0.55*D25=262; no nucleus
; ###########################################################################
; D25=151; 0.55*D25=82.9; crummy; no nucleus
; jm08jan15nyu - crappy!
;   ispec, 'holmberg_ix_drift_020.fits', /deredden, aperture=drift20, $
;     /meanprofile, tracespec=1, traceorder=1, s1, $
;     datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; ###########################################################################
; D25=791; 0.55*D25=435; no nucleus
    ispec, 'ic_2574_drift_020.fits', /deredden, aperture=drift20, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; ###########################################################################
; IC4710 - no data; D25=218; 0.55*D25=120
; ###########################################################################
; jm08jan15nyu - crappy!
; crummy
;   ispec, 'm81_dwa_drift_020.fits', /deredden, aperture=drift20, $
;     /meanprofile, tracespec=1, traceorder=1, s1, $
;     datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; D25=78; 0.55*D25=42.9; foreground star subtracted
;   ispec, 'm81_dwa_drift_055.fits', /deredden, refrow=52, aperture=42.9, $
;     /meanprofile, tracespec=1, traceorder=1, s1, $
;     /starfit, star_refrow=72, debug=debug, $
;     datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; ###########################################################################
    ispec, 'm81_dwb_drift_020.fits', /deredden, aperture=drift20, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; D25=52.3; 0.55*D25=28.8
    ispec, 'm81_dwb_drift_055.fits', /deredden, aperture=28.8, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; ###########################################################################
    ispec, 'mrk_0033_drift_020.fits', /deredden, aperture=drift20, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; D25=60.0; 0.55*D25=33.0
    ispec, 'mrk_0033_drift_055.fits', /deredden, aperture=33.0, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
    ispec, 'mrk_0033_nuclear.fits', /deredden, aperture=nucap, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; ###########################################################################
; extended spatial profile    
    ispec, 'ngc_0024_drift_020.fits', /deredden, aperture=drift20, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; D25=345; 0.55*D25=190
    ispec, 'ngc_0024_drift_055.fits', /deredden, aperture=190.0, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
    ispec, 'ngc_0024_nuclear.fits', /deredden, aperture=nucap, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; ###########################################################################
; great    
    ispec, 'ngc_0337_drift_020.fits', /deredden, aperture=drift20, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; D25=173; 0.55*D25=95.2
    ispec, 'ngc_0337_drift_055.fits', /deredden, aperture=95.2, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
    ispec, 'ngc_0337_nuclear.fits', /deredden, aperture=nucap, $
      /meanprofile, tracespec=1, traceorder=1, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; ###########################################################################
; problems with the telluric correction    
    ispec, 'ngc_0584_drift_020.fits', /deredden, aperture=drift20, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; D25=250; 0.55*D25=138
    ispec, 'ngc_0584_drift_055.fits', /deredden, aperture=138.0, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
    ispec, 'ngc_0584_nuclear.fits', /deredden, aperture=nucap, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; ###########################################################################
    ispec, 'ngc_0628_drift_020.fits', /deredden, aperture=drift20, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
    ispec, 'ngc_0628_nuclear.fits', /deredden, aperture=nucap, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; D25=628; 0.55*D25=346; foreground star subtracted
    ispec, 'ngc_0628_drift_055.fits', /deredden, aperture=346.0, $
       /meanprofile, tracespec=1, inputspec=inputspec, datapath=datapath, outpath=outpath, $
      /starfit, star_refrow=60, center_order=1, sigma_order=1, debug=debug, $
      noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; ###########################################################################
    ispec, 'ngc_0855_drift_020.fits', /deredden, aperture=drift20, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; D25=158; 0.55*D25=86.8
    ispec, 'ngc_0855_drift_055.fits', /deredden, aperture=86.8, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
    ispec, 'ngc_0855_nuclear.fits', /deredden, aperture=nucap, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; ###########################################################################
    ispec, 'ngc_0925_drift_020.fits', /deredden, refrow=65, aperture=drift20, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
    ispec, 'ngc_0925_nuclear.fits', /deredden, aperture=nucap, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; D25=628; 0.55*D25=346; foreground stars subtracted; unable to
; accurately subtract the star near row 175
    ispec, 'ngc_0925_drift_055.fits', /deredden, refrow=112, aperture=346.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      /starfit, star_refrow=[80,200], sigma_order=1, debug=debug, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; ###########################################################################
; difficult traces on the 20" and nuclear spectra due to a double
; nucleus???     
    ispec, 'ngc_1097_drift_020.fits', /deredden, aperture=drift20, $
      /meanprofile, tracespec=1, traceorder=1, trace_mincol=600, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; D25=560; 0.55*D25=308
    ispec, 'ngc_1097_drift_056.fits', /deredden, aperture=308.0, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
    ispec, 'ngc_1097_nuclear_1.fits', /deredden, refrow=107, aperture=nucap, $
      /meanprofile, tracespec=1, traceorder=1, trace_mincol=200, s1, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
    ispec, 'ngc_1097_nuclear_2.fits', /deredden, refrow=107, aperture=nucap, $
      /meanprofile, tracespec=1, traceorder=1, trace_mincol=500, s2, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; ###########################################################################
    ispec, 'ngc_1266_drift_020.fits', /deredden, aperture=drift20, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; D25=92.9; 0.55*D25=51.1
    ispec, 'ngc_1266_drift_055.fits', /deredden, refrow=65, aperture=51.1, $
      /meanprofile, fluxrange=[-0.2,2]*1E-15, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
    ispec, 'ngc_1266_nuclear.fits', /deredden, aperture=nucap, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; ###########################################################################
; repeat observations; virtually *no* signal in the second spectrum
; due to clouds, so just use the first spectrum
    ispec, 'ngc_1291_drift_1_020.fits', /deredden, refrow=110, aperture=drift20, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric, $
      outname='ngc_1291_drift_020'
    ispec, 'ngc_1291_drift_2_020.fits', /deredden, refrow=110, aperture=drift20, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=0, /nobadmask, /noskymask, telluric=telluric
; D25=586; 0.55*D25=323; slit length = 308"; foreground star
; subtracted; only extract 293" (=0.5*D25)
    ispec, 'ngc_1291_drift_056.fits', /deredden, aperture=293.0, $
      /meanprofile, tracespec=1, traceorder=1, $
      /starfit, star_refrow=28, debug=debug, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
    ispec, 'ngc_1291_nuclear.fits', /deredden, refrow=110, aperture=nucap, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; ###########################################################################
; great    
    ispec, 'ngc_1316_drift_020.fits', /deredden, aperture=drift20, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; D25=721; 0.55*D25=397; slit length = 308"; only extract 288"
; (=0.4*D25) 
    ispec, 'ngc_1316_drift_056.fits', /deredden, aperture=288.0, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
    ispec, 'ngc_1316_nuclear.fits', /deredden, aperture=nucap, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; ###########################################################################
; good    
    ispec, 'ngc_1377_drift_020.fits', /deredden, aperture=drift20, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; D25=107; 0.55*D25=58.7
    ispec, 'ngc_1377_drift_055.fits', /deredden, aperture=58.7, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
    ispec, 'ngc_1377_nuclear.fits', /deredden, aperture=nucap, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; ###########################################################################
; great    
    ispec, 'ngc_1404_drift_020.fits', /deredden, aperture=drift20, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; D25=199; 0.55*D25=109
    ispec, 'ngc_1404_drift_056.fits', /deredden, aperture=109.0, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
    ispec, 'ngc_1404_nuclear.fits', /deredden, aperture=nucap, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; ###########################################################################
; great    
    ispec, 'ngc_1482_drift_020.fits', /deredden, aperture=drift20, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
    ispec, 'ngc_1482_nuclear.fits', /deredden, aperture=nucap, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; D25=147; 0.55*D25=81.0; repeat observations
    ispec, 'ngc_1482_drift_055.fits', /deredden, aperture=81.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
    ispec, 'ngc_1482_drift_1_056.fits', /deredden, aperture=81.0, $
      /meanprofile, tracespec=1, traceorder=1, s2, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
    ispec, 'ngc_1482_drift_2_056.fits', /deredden, aperture=81.0, $
      /meanprofile, tracespec=1, traceorder=1, s3, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; since it appears that we got a nuclear spectrum (above) do not
; extract these east and west "knot" spectra    
    ispec, 'ngc_1482_e_knot.fits', /deredden, aperture=nucap, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=0, /nobadmask, /noskymask, telluric=telluric
    ispec, 'ngc_1482_w_knot.fits', /deredden, aperture=nucap, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=0, /nobadmask, /noskymask, telluric=telluric
; ###########################################################################
; fantastic    
    ispec, 'ngc_1512_drift_020.fits', /deredden, aperture=drift20, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; D25=535; 0.55*D25=294
    ispec, 'ngc_1512_drift_056.fits', /deredden, aperture=294.0, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
    ispec, 'ngc_1512_nuclear.fits', /deredden, aperture=nucap, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; ###########################################################################
; fantastic    
    ispec, 'ngc_1566_drift_020.fits', /deredden, aperture=drift20, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; D25=499; 0.55*D25=275
    ispec, 'ngc_1566_drift_056.fits', /deredden, aperture=275.0, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
    ispec, 'ngc_1566_nuclear.fits', /deredden, aperture=nucap, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; ###########################################################################
; very blue nucleus!!    
    ispec, 'ngc_1705_drift_020.fits', /deredden, aperture=drift20, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; D25=114; 0.55*D25=62.9
    ispec, 'ngc_1705_drift_056.fits', /deredden, aperture=62.9, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
    ispec, 'ngc_1705_nuclear.fits', /deredden, aperture=nucap, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; ###########################################################################
; D25=1313; 0.55*D25=722; star or not?
    ispec, 'ngc_2403_drift_056.fits', /deredden, refrow=300, aperture=657.0, $; 722.0, $
      /meanprofile, tracespec=1, traceorder=1, $
      /starfit, star_refrow=[289,215,414,523], frac_rows=0.03, ngauss_terms=6L, debug=0, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; extended spatial profile
    ispec, 'ngc_2403_drift_020.fits', /deredden, aperture=drift20, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
;   ispec, 'ngc_2403_nuclear.fits', /deredden, aperture=nucap, $
;     /meanprofile, tracespec=1, traceorder=1, $
;     datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; ###########################################################################
; very good; can also get NGC2799
    ispec, 'ngc_2798_drift_020.fits', /deredden, aperture=drift20, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; D25=154; 0.55*D25=84.8
    ispec, 'ngc_2798_drift_055.fits', /deredden, aperture=84.8, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
    ispec, 'ngc_2798_nuclear.fits', /deredden, aperture=nucap, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; ###########################################################################
; D25=488; 0.55*D25=268
    ispec, 'ngc_2841_drift_056.fits', /deredden, refrow=120, aperture=268.0, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; great    
    ispec, 'ngc_2841_drift_020.fits', /deredden, aperture=drift20, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
    ispec, 'ngc_2841_nuclear.fits', /deredden, aperture=nucap, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; ###########################################################################
; great    
    ispec, 'ngc_2915_drift_020.fits', /deredden, aperture=drift20, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; D25=114; 0.55*D25=62.9
    ispec, 'ngc_2915_drift_056.fits', /deredden, aperture=62.9, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
    ispec, 'ngc_2915_nuclear.fits', /deredden, aperture=nucap, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; ###########################################################################
; D25=353; 0.55*D25=194
    ispec, 'ngc_2976_drift_056.fits', /deredden, aperture=194.0, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; extended spatial profile    
    ispec, 'ngc_2976_drift_020.fits', /deredden, aperture=drift20, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
    ispec, 'ngc_2976_nuclear.fits', /deredden, aperture=nucap, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; ###########################################################################
; D25=1615; 0.55*D25=888; foreground star subtracted
    ispec, 'ngc_3031_drift_056.fits', /deredden, aperture=800.0, $
      /meanprofile, tracespec=1, traceorder=1, $
      /starfit, star_refrow=250, debug=debug, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; great 
    ispec, 'ngc_3031_drift_020.fits', /deredden, aperture=drift20, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
    ispec, 'ngc_3031_nuclear.fits', /deredden, aperture=nucap, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; ###########################################################################
; D25=673; 0.55*D25=370; offset from nucleus, so use a dedicated
; extraction aperture
; jm08jan15nyu - the radial strip spectrum is not very useful
;   ispec, 'ngc_3034_drift_056.fits', /deredden, aperture=120.0, $
;     /meanprofile, tracespec=1, traceorder=1, $
;     datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
    ispec, 'ngc_3034_drift_020.fits', /deredden, aperture=drift20, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
    ispec, 'ngc_3034_nuclear.fits', /deredden, aperture=nucap, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; ###########################################################################
; great    
    ispec, 'ngc_3049_drift_020.fits', /deredden, aperture=drift20, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; D25=131; 0.55*D25=72.2
    ispec, 'ngc_3049_drift_055.fits', /deredden, aperture=72.2, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
    ispec, 'ngc_3049_nuclear.fits', /deredden, aperture=nucap, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; ###########################################################################
; good
    ispec, 'ngc_3184_drift_020.fits', /deredden, aperture=drift20, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; D25=445; 0.55*D25=245
    ispec, 'ngc_3184_drift_055.fits', /deredden, aperture=245.0, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
    ispec, 'ngc_3184_nuclear.fits', /deredden, aperture=nucap, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; ###########################################################################
; great
    ispec, 'ngc_3190_drift_020.fits', /deredden, aperture=drift20, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; D25=262; 0.55*D25=144
    ispec, 'ngc_3190_drift_055.fits', /deredden, aperture=144.0, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
    ispec, 'ngc_3190_nuclear.fits', /deredden, aperture=nucap, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; ###########################################################################
; great
    ispec, 'ngc_3198_drift_020.fits', /deredden, aperture=drift20, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; D25=511; 0.55*D25=281; only extract 180" (=0.35*D25)
    ispec, 'ngc_3198_drift_055.fits', /deredden, aperture=180.0, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
    ispec, 'ngc_3198_nuclear.fits', /deredden, aperture=nucap, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; ###########################################################################
; great; small residual CR interpolation problems in the blue
    ispec, 'ngc_3265_drift_020.fits', /deredden, refrow=60, aperture=drift20, $
      /meanprofile, tracespec=1, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; D25=77.3; 0.55*D25=42.5; foreground star subtracted
    ispec, 'ngc_3265_drift_055.fits', /deredden, aperture=42.5, $
      /meanprofile, tracespec=1, traceorder=1, $
      /starfit, star_refrow=42, debug=debug, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
    ispec, 'ngc_3265_nuclear.fits', /deredden, aperture=nucap, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; ###########################################################################
; repeat observations; the first scan has higher S/N
    ispec, 'ngc_3351_drift_1_020.fits', /deredden, aperture=drift20, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
    ispec, 'ngc_3351_drift_2_020.fits', /deredden, aperture=drift20, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; D25=445; 0.55*D25=245; the KPNO spectrum has higher S/N    
    ispec, 'ngc_3351_drift_056.fits', /deredden, aperture=245.0, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
    ispec, 'ngc_3351_drift_055.fits', /deredden, aperture=245.0, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; double nucleus??; the first spectrum was taken on 2001-Dec at an
; airmass of 1.54 and a position angle of 110; the parallactic angle
; was 27 degrees at the time, so there is lots of blue light lost.  so
; choose the KPNO (second) spectrum
    ispec, 'ngc_3351_nuclear_1.fits', /deredden, aperture=nucap, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=0, /nobadmask, /noskymask, telluric=telluric
    ispec, 'ngc_3351_nuclear_2.fits', /deredden, aperture=nucap, $
      /meanprofile, tracespec=1, traceorder=1, s2, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric, $
      outname='ngc_3351_nuclear'
; ###########################################################################
; repeat observations; the first scan has higher S/N
    ispec, 'ngc_3521_drift_1_020.fits', /deredden, aperture=drift20, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
    ispec, 'ngc_3521_drift_2_020.fits', /deredden, aperture=drift20, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; D25=658; 0.55*D25=362; only extract 263" (=0.4*D25); slit
; length=307" 
    ispec, 'ngc_3521_drift_056.fits', /deredden, aperture=263.0, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
    ispec, 'ngc_3521_nuclear.fits', /deredden, aperture=nucap, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; ###########################################################################
; great    
    ispec, 'ngc_3621_drift_020.fits', /deredden, refrow=110, aperture=drift20, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; D25=738; 0.55*D25=406; foreground stars subtracted; asymmetric
; spatial profile; only extract 230" (=0.31*D25)
    ispec, 'ngc_3621_drift_056.fits', /deredden, refrow=110, aperture=230.0, $
      /meanprofile, tracespec=1, traceorder=1, $
      /starfit, star_refrow=[201,191], trace_searchbox=5, debug=debug, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
    ispec, 'ngc_3621_nuclear.fits', /deredden, aperture=nucap, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; ###########################################################################
; great    
    ispec, 'ngc_3627_drift_020.fits', /deredden, aperture=drift20, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; D25=547; 0.55*D25=301; only extract 200" (=0.37*D25)
    ispec, 'ngc_3627_drift_055.fits', /deredden, aperture=200.0, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
    ispec, 'ngc_3627_nuclear.fits', /deredden, aperture=nucap, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; ###########################################################################
; great    
    ispec, 'ngc_3773_drift_020.fits', /deredden, aperture=drift20, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; D25=70.5; 0.55*D25=38.8
    ispec, 'ngc_3773_drift_055.fits', /deredden, aperture=38.8, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
    ispec, 'ngc_3773_nuclear.fits', /deredden, aperture=nucap, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; ###########################################################################
; great    
    ispec, 'ngc_3938_drift_020.fits', /deredden, aperture=drift20, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; D25=322; 0.55*D25=177; pretty crummy
    ispec, 'ngc_3938_drift_056.fits', /deredden, aperture=177.0, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
    ispec, 'ngc_3938_nuclear.fits', /deredden, aperture=nucap, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; ###########################################################################
; great    
    ispec, 'ngc_4125_drift_020.fits', /deredden, aperture=drift20, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; D25=345; 0.55*D25=190
    ispec, 'ngc_4125_drift_056.fits', /deredden, aperture=190.0, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
    ispec, 'ngc_4125_nuclear.fits', /deredden, aperture=nucap, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; ###########################################################################
    ispec, 'ngc_4236_drift_020.fits', /deredden, aperture=drift20, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; D25=1312; 0.55*D25=722; only extract 656" (=0.5*D25); do not trace;
; foreground star subtracted; crummy; remove! crappy!
;   ispec, 'ngc_4236_drift_056.fits', /deredden, refrow=209, aperture=656.0, $
;     /meanprofile, tracespec=0, traceorder=1, $
;     /starfit, star_refrow=180, debug=debug, $
;     datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; ###########################################################################
; great    
    ispec, 'ngc_4254_drift_020.fits', /deredden, aperture=drift20, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; D25=322; 0.55*D25=177
    ispec, 'ngc_4254_drift_055.fits', /deredden, aperture=177.0, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
    ispec, 'ngc_4254_nuclear.fits', /deredden, aperture=nucap, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; ###########################################################################
; great    
    ispec, 'ngc_4321_drift_020.fits', /deredden, aperture=drift20, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; D25=445; 0.55*D25=245
    ispec, 'ngc_4321_drift_056.fits', /deredden, aperture=245.0, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
    ispec, 'ngc_4321_nuclear.fits', /deredden, aperture=nucap, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; ###########################################################################
; great    
    ispec, 'ngc_4450_drift_020.fits', /deredden, aperture=drift20, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; D25=315; 0.55*D25=173
    ispec, 'ngc_4450_drift_055.fits', /deredden, aperture=173.0, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
    ispec, 'ngc_4450_nuclear.fits', /deredden, aperture=nucap, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; ###########################################################################
; great    
    ispec, 'ngc_4536_drift_020.fits', /deredden, aperture=drift20, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; D25=455; 0.55*D25=250
    ispec, 'ngc_4536_drift_055.fits', /deredden, aperture=250.0, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
    ispec, 'ngc_4536_nuclear.fits', /deredden, aperture=nucap, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; ###########################################################################
; great    
    ispec, 'ngc_4552_drift_020.fits', /deredden, aperture=drift20, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; D25=308; 0.55*D25=169
    ispec, 'ngc_4552_drift_055.fits', /deredden, aperture=169.0, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
    ispec, 'ngc_4552_nuclear.fits', /deredden, aperture=nucap, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; ###########################################################################
; repeat observations
    ispec, 'ngc_4559_drift_020.fits', /deredden, refrow=60, aperture=drift20, $
      /meanprofile, tracespec=1, traceorder=1, fluxrange=[0,150]*1E-17, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; D25=643; 0.55*D25=354; foreground star subtracted
    ispec, 'ngc_4559_drift_055.fits', /deredden, refrow=115, aperture=354.0, $
      /meanprofile, tracespec=1, traceorder=1, fluxrange=[-10,250]*1E-17, $
      /starfit, star_refrow=35, debug=debug, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
    ispec, 'ngc_4559_nuclear_1.fits', /deredden, refrow=60, aperture=nucap, $
      /meanprofile, tracespec=1, traceorder=1, fluxrange=[0,40]*1E-17, s1, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
    ispec, 'ngc_4559_nuclear_2.fits', /deredden, refrow=60, aperture=nucap, $
      /meanprofile, tracespec=1, traceorder=1, fluxrange=[0,40]*1E-17, s2, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; ###########################################################################
; great    
    ispec, 'ngc_4569_drift_020.fits', /deredden, aperture=drift20, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; D25=573; 0.55*D25=315; only extracted 172" (=0.3*D25)
    ispec, 'ngc_4569_drift_055.fits', /deredden, aperture=172.0, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
    ispec, 'ngc_4569_nuclear.fits', /deredden, aperture=nucap, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; ###########################################################################
; repeat observations; the second scan has higher S/N
    ispec, 'ngc_4579_drift_1_020.fits', /deredden, aperture=drift20, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
    ispec, 'ngc_4579_drift_2_020.fits', /deredden, aperture=drift20, $
      /meanprofile, tracespec=1, traceorder=1, s2, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; D25=353; 0.55*D25=194
    ispec, 'ngc_4579_drift_055.fits', /deredden, aperture=194.0, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
    ispec, 'ngc_4579_nuclear.fits', /deredden, aperture=nucap, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; ###########################################################################
; repeat observations; the second scan has higher S/N
    ispec, 'ngc_4594_drift_1_020.fits', /deredden, aperture=drift20, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
    ispec, 'ngc_4594_drift_2_020.fits', /deredden, aperture=drift20, $
      /meanprofile, tracespec=1, traceorder=1, s2, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; D25=523; 0.55*D25=287; Sombrero
    ispec, 'ngc_4594_drift_055.fits', /deredden, aperture=287.0, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; the second scan has higher S/N
    ispec, 'ngc_4594_nuclear_1.fits', /deredden, aperture=nucap, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
    ispec, 'ngc_4594_nuclear_2.fits', /deredden, aperture=nucap, $
      /meanprofile, tracespec=1, traceorder=1, s2, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; ###########################################################################
; good    
    ispec, 'ngc_4625_drift_020.fits', /deredden, aperture=drift20, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; D25=131; 0.55*D25=72.2
    ispec, 'ngc_4625_drift_056.fits', /deredden, aperture=72.2, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
    ispec, 'ngc_4625_nuclear.fits', /deredden, aperture=nucap, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; ###########################################################################
; no nucleus
; D25=930; 0.55*D25=512; do not bother tracing
    ispec, 'ngc_4631_drift_056.fits', /deredden, refrow=180, aperture=512.0, $
      /meanprofile, tracespec=0, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; extended spatial profile, centered on the HII Region
    ispec, 'ngc_4631_drift_020.fits', /deredden, refrow=40.0, aperture=drift20, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; ###########################################################################
; great    
    ispec, 'ngc_4725_drift_020.fits', /deredden, aperture=drift20, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; D25=929; 0.55*D25=511; foreground star subtracted; only extract 292"
; (=0.31*D25) 
    ispec, 'ngc_4725_drift_055.fits', /deredden, aperture=292.0, $
      /meanprofile, tracespec=1, traceorder=1, $
      /starfit, star_refrow=150, frac_rows=0.06, sigma_order=1, debug=debug, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
    ispec, 'ngc_4725_nuclear.fits', /deredden, aperture=nucap, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; ###########################################################################
; great    
    ispec, 'ngc_4736_drift_020.fits', /deredden, aperture=drift20, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; D25=673; 0.55*D25=370
    ispec, 'ngc_4736_drift_056.fits', /deredden, aperture=370.0, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
    ispec, 'ngc_4736_nuclear.fits', /deredden, aperture=nucap, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; ###########################################################################
; great    
    ispec, 'ngc_4826_drift_020.fits', /deredden, aperture=drift20, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; D25=600; 0.55*D25=330
    ispec, 'ngc_4826_drift_055.fits', /deredden, aperture=330.0, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
    ispec, 'ngc_4826_nuclear.fits', /deredden, aperture=nucap, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; ###########################################################################
; D25=643; 0.55*D25=354; foreground star subtracted
    ispec, 'ngc_5033_drift_056.fits', /deredden, aperture=354.0, $
      /meanprofile, tracespec=1, traceorder=1, $
      /starfit, star_refrow=105, debug=debug, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; great    
    ispec, 'ngc_5033_drift_020.fits', /deredden, aperture=drift20, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
    ispec, 'ngc_5033_nuclear.fits', /deredden, aperture=nucap, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; ###########################################################################
; great    
    ispec, 'ngc_5055_drift_020.fits', /deredden, aperture=drift20, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; D25=755; 0.55*D25=416
    ispec, 'ngc_5055_drift_056.fits', /deredden, aperture=416.0, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
    ispec, 'ngc_5055_nuclear.fits', /deredden, aperture=nucap, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; ###########################################################################
; possible poor CR interpolation near ~5900 Angstroms (see 02feb - cra.3856)
    ispec, 'ngc_5194_drift_020.fits', /deredden, aperture=drift20, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; D25=673; 0.55*D25=370
    ispec, 'ngc_5194_drift_056.fits', /deredden, aperture=370.0, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; great    
    ispec, 'ngc_5194_nuclear.fits', /deredden, aperture=nucap, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; ###########################################################################
; great; dusty!
    ispec, 'ngc_5195_drift_020.fits', /deredden, aperture=drift20, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; D25=345; 0.55*D25=190
    ispec, 'ngc_5195_drift_056.fits', /deredden, aperture=190.0, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
    ispec, 'ngc_5195_nuclear.fits', /deredden, aperture=nucap, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; ###########################################################################
; Tol89=NGC5398; okay 
    ispec, 'tololo_89_drift_020.fits', /deredden, aperture=drift20, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; D25=169; 0.55*D25=93.0; stellar contamination
    ispec, 'tololo_89_drift_055.fits', /deredden, refrow=65, aperture=93.0, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; ###########################################################################
; NGC5408 - no data; D25=97.3; 0.55*D25=53.5
; ###########################################################################
; D25=287; 0.55*D25=158
    ispec, 'ngc_5474_drift_056.fits', /deredden, aperture=158.0, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; good
    ispec, 'ngc_5474_drift_020.fits', /deredden, aperture=drift20, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; ###########################################################################
; repeat observations; D25=165; 0.55*D25=90.0; great; the first
; (PA=113) spectrum has higher S/N than the second one (PA=90)
    ispec, 'ngc_5713_drift_055.fits', /deredden, aperture=90.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
    ispec, 'ngc_5713_drift_056.fits', /deredden, aperture=90.0, $
      /meanprofile, tracespec=1, traceorder=1, s2, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
    ispec, 'ngc_5713_drift_020.fits', /deredden, aperture=drift20, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
    ispec, 'ngc_5713_nuclear.fits', /deredden, aperture=nucap, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; ###########################################################################
; great    
    ispec, 'ngc_5866_drift_020.fits', /deredden, aperture=drift20, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; D25=281; 0.55*D25=154
    ispec, 'ngc_5866_drift_055.fits', /deredden, aperture=154.0, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
    ispec, 'ngc_5866_nuclear.fits', /deredden, aperture=nucap, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; ###########################################################################
; NGC6822 - no data; D25=929; 0.55*D25=511
; ###########################################################################
; D25=689; 0.55*D25=379; foreground stars subtracted; stellar
; contamination 
    ispec, 'ngc_6946_drift_056.fits', /deredden, refrow=220, aperture=379.0, $
      /meanprofile, tracespec=1, traceorder=1, $
      /starfit, star_refrow=[303,427], frac_row=0.05, debug=0, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; dusty!  
    ispec, 'ngc_6946_drift_020.fits', /deredden, aperture=drift20, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
    ispec, 'ngc_6946_nuclear.fits', /deredden, aperture=nucap, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; ###########################################################################
; great    
    ispec, 'ngc_7331_drift_020.fits', /deredden, aperture=drift20, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; D25=628; 0.55*D25=346; only extract 200" (=0.31*D25) 
    ispec, 'ngc_7331_drift_055.fits', /deredden, aperture=200.0, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
    ispec, 'ngc_7331_nuclear.fits', /deredden, aperture=nucap, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; ###########################################################################
; NGC7552 - no data; D25=203; 0.55*D25=112 
; ###########################################################################
; great    
    ispec, 'ngc_7793_drift_020.fits', /deredden, aperture=drift20, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; D25=560; 0.55*D25=308
    ispec, 'ngc_7793_drift_056.fits', /deredden, aperture=308.0, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; very blue!  post-burst nucleus!
    ispec, 'ngc_7793_nuclear.fits', /deredden, aperture=nucap, $
      /meanprofile, tracespec=1, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, telluric=telluric
; ###########################################################################

; combine repeaters

    splog, 'Combining repeaters.'
    sings_combine_repeaters, wfits=wfits
    
    splog, 'Total time to extract all SINGS spectra (minutes):', (systime(1)-t0)/60.0

; sings_header_redshifts, /update    
; sings_info, info, /write
;
; sings_specfit, /nuclear, /drift20, /drift56

; sings_sigspecrepair, /nuclear, /repair
; sings_sigspecrepair, /drift20, /repair
; sings_sigspecrepair, /drift56, /repair
    
; sings_skyrepair, /nuclear, /repair
; sings_skyrepair, /drift20, /repair
; sings_skyrepair, /drift56, /repair
    
return
end
