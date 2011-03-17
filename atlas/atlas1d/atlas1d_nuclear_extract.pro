pro atlas1d_nuclear_extract, wfits=wfits, noplot=noplot
; J. Moustakas 2003 Dec 29 13:26

; Excluded nuclear galaxies (2005-Jul-30):
; ----------------------------------------

; NGC1614 - the spectrum is simply wrong, even though the reductions
;           look fine; so exclude from the atlas (jm05jul24uofa)

; NGC1560 - poorly-defined nucleus, so exclude (jm05jul24uofa)    

; IRAS15250+3690 - H-alpha is redshifted beyond the red wavelength
;                  cutoff     
    
; Repeat observations:
; --------------------
; NGC0877
; NGC1058
; NGC1275
; NGC2541
; NGC2599
; NGC3049
; NGC4559
; NGC4594
; NGC6240

; Interacting pairs:
; -----------------
; NGC5929/NGC5930 - separate (nuclear)

; Difficult sky-subtraction:
; --------------------------
; NGC1003, NGC1421, NGC2276, NGC2500, NGC2541, NGC3079, NGC3628,
; NGC3917, NGC4254, NGC4559, UGCA114

    datapath = atlas_path(/atlas2d)
    outpath = atlas_path(/atlas1d)
    repeatpath = atlas_path(/atlas1d)+'repeaters/'

    nucap = 2.5 ; nuclear aperture [arcsec]

    t0 = systime(1)

    ispec, 'ic_0750_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; Turner    
    ispec, 'ic_0883_nuclear.fits', /deredden, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits
    
    ispec, 'ic_1076_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; IRAS15250+3609 has been excluded from the atlas (02apr) because
; H-alpha is redshifted beyond the spectral coverage
    ispec, 'iras_15250+3609_nuclear.fits', /deredden, aperture=nucap, $
      /meanprofile, tracespec=1, traceorder=2, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=0, /noskymask, telluric=telluric
    
    ispec, 'mrk_0055_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'mrk_0315_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'mrk_0685_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'mrk_0848_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'mrk_1460_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'mrk_1490_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_0024_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_0337_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_0584_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_0615_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_0628_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_0694_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_0695_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_0855_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; ###########################################################################
    ispec, 'ngc_0877_nuclear_1.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, s1, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_0877_nuclear_2.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, s2, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric
; ###########################################################################

    ispec, 'ngc_0925_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_0958_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_0972_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_0976_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_1003_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_1036_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; ###########################################################################
; comparable S/N = 21 in both observations
    ispec, 'ngc_1058_nuclear_1.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_1058_nuclear_2.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric
; ###########################################################################

    ispec, 'ngc_1068_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_1084_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_1087_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_1134_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_1143_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_1144_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; poorly-defined nucleus; center on an HII region
    ispec, 'ngc_1156_nuclear.fits', /deredden, refrow=70, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_1266_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; ###########################################################################
; comparable S/N = 38 in both observations
    ispec, 'ngc_1275_nuclear_1.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_1275_nuclear_2.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric
; ###########################################################################

    ispec, 'ngc_1345_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_1357_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_1359_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_1377_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_1385_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_1421_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; poorly-defined nucleus; exclude
    ispec, 'ngc_1560_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, fluxrange=[0,15]*1E-17, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=0, /noskymask, telluric=telluric

; wrong, but I don't know why! (jm05jul24uofa)    
    ispec, 'ngc_1614_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=0, /noskymask, telluric=telluric

    ispec, 'ngc_1800_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_1832_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_2090_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_2139_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_2146_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_2276_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_2337_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_2388_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_2415_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_2500_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; ###########################################################################
; comparably low S/N = 7 in both observations    
    ispec, 'ngc_2541_nuclear_1.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_2541_nuclear_2.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluricb
; ###########################################################################

; ###########################################################################
; S/N = 45, 31, and 47
    ispec, 'ngc_2599_nuclear_1.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_2599_nuclear_2.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_2599_nuclear_3.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric
; ###########################################################################

    ispec, 'ngc_2775_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_2782_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_2798_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_2841_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_2893_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_2903_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_3031_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; ###########################################################################
; comparable S/N = 18 in both observations
    ispec, 'ngc_3049_nuclear_1.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_3049_nuclear_2.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric
; ###########################################################################

    ispec, 'ngc_3079_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_3184_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_3190_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_3198_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_3239_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_3265_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_3274_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; possible CR sitting on H-alpha    
    ispec, 'ngc_3310_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_3344_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_3351_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_3353_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_3365_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_3395_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_3396_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_3442_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_3504_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_3521_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_3600_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_3627_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_3628_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; extremely dusty!!!!
    ispec, 'ngc_3718_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_3726_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_3729_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; poorly-defined nucleus
    ispec, 'ngc_3769_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_3773_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_3782_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_3870_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_3896_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; poorly defined nucleus
    ispec, 'ngc_3906_nuclear.fits', /deredden, refrow=70, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=5, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_3913_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_3917_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; Turner    
    ispec, 'ngc_3921_nuclear.fits', /deredden, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits

    ispec, 'ngc_3928_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_3938_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_3949_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_3953_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_3972_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_3991_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_3994_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_3995_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_3998_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_4010_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_4051_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_4088_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_4100_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_4102_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_4111_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_4125_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_4136_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_4138_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_4144_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_4150_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_4163_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_4190_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_4217_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_4218_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_4220_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_4254_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_4258_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_4303_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_4321_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_4384_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_4389_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_4414_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; mistakenly called NGC 4550 previously! (jm04oct12uofa)
    ispec, 'ngc_4450_nuclear.fits', /deredden, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_4500_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_4536_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_4552_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; ###########################################################################
; the second observation is at higher airmass and is not precisely at
; the parallactic angle; S/N = 16 and 15
    ispec, 'ngc_4559_nuclear_1.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, fluxrange=[0,40]*1E-17, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_4559_nuclear_2.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, fluxrange=[0,40]*1E-17, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric
; ###########################################################################

    ispec, 'ngc_4569_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_4579_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; ###########################################################################
; although the second exposure has the larger S/N and both objects
; were observed at high airmass, the first exposure was taken at
; precisely the parallactic angle, while the second exposure was ~7
; degrees off; choose the first exposure

    ispec, 'ngc_4594_nuclear_1.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_4594_nuclear_2.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric
; ###########################################################################

    ispec, 'ngc_4625_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_4670_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; Turner    
    ispec, 'ngc_4676a_nuclear.fits', /deredden, refrow=100.0, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=1, sbox=10, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits
    
; Turner    
    ispec, 'ngc_4676b_nuclear.fits', /deredden, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits

    ispec, 'ngc_4725_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_4736_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_4826_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_4900_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_5014_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_5033_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_5055_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_5104_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_5144_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_5194_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_5195_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_5430_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_5548_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_5653_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_5676_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_5713_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_5866_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; ###########################################################################
; interacting pair with NGC5930;  44-degree position angle; NGC5930 is
; NE of 5929
    ispec, 'ngc_5929_nuclear.fits', /deredden, refrow=80, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='ngc_5929_nuclear', galaxy='NGC5929'
    
    ispec, 'ngc_5929_nuclear.fits', /deredden, refrow=50, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='ngc_5930_nuclear', galaxy='NGC5930'
; ###########################################################################
    
    ispec, 'ngc_5936_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_5953_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_5992_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_5996_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_6052_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_6090_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_6181_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_6207_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_6217_nuclear.fits', /deredden, refrow=80, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; ###########################################################################
; comparable S/N = 13 on both observations    
    ispec, 'ngc_6240_nuclear_1.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_6240_nuclear_2.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric
; ###########################################################################

    ispec, 'ngc_6285_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_6286_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_6701_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_6946_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_7137_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_7244_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_7331_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_7448_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_7465_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_7468_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_7518_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_7678_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_7782_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ugc_00685_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ugc_01561_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ugc_02238_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ugc_02982_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ugc_03426_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ugc_05151_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ugc_05720_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ugc_06103_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ugc_06917_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ugc_07354_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ugc_09081_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ugca_073_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ugca_090_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ugca_114_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ugca_116_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ugca_410_nuclear.fits', /deredden, refrow=62, aperture=nucap, $
      /meanprofile, /tracespec, traceorder=2, sbox=10, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    splog, 'Total time to extract all nuclear spectra (minutes):', (systime(1)-t0)/60.0

; write redshifts

;   atlas1d_redshifts, /update

return
end
