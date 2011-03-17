pro turner_mergers_extract, wfits=wfits, noplot=noplot

;   datapath = atlas_path(/spec2dturner)+'spec2d/'
    datapath = atlas_path(/spec2dturner)+'spec2d/rectified/'
    outpath = atlas_path(/spec2dturner)+'spec1d/'
    repeatpath = atlas_path(/spec2dturner)+'spec1d/repeaters/'

; ##################################################    
; NUCLEAR
; ##################################################    
    
    ispec, 'ic_0883_nuclear.fits', /deredden, aperture=2.5, $
      /meanprofile, /tracespec, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits
    
    ispec, 'ngc_3921_nuclear.fits', /deredden, aperture=2.5, $
      /meanprofile, /tracespec, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits

    ispec, 'ngc_4676a_nuclear.fits', /deredden, refrow=100.0, aperture=2.5, $
      /meanprofile, /tracespec, traceorder=1, sbox=10, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits
    
    ispec, 'ngc_4676b_nuclear.fits', /deredden, aperture=2.5, $
      /meanprofile, /tracespec, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits
    
; ##################################################    
; INTEGRATED
; ##################################################    
    
; ARP230; the object on the shoulder is a residual cosmic ray; the
; trace fails, so do not trace; Turner
    ispec, 'ic_0051_drift_089.fits', /deredden, aperture=75.0, $
      /meanprofile, tracespec=0, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits

; ARP228 = UGC01267
    ispec, 'ic_0162_drift_040.fits', /deredden, aperture=100.0, $
      /meanprofile, /tracespec, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits

; ###########################################################################    
; repeat observations
    ispec, 'ic_0883_drift_1_020.fits', /deredden, aperture=90.0, $
      /meanprofile, /tracespec, traceorder=1, s1, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits

    ispec, 'ic_0883_drift_2_020.fits', /deredden, aperture=90.0, $
      /meanprofile, /tracespec, traceorder=1, s2, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits
; ###########################################################################    

    ispec, 'ngc_0474_drift_060.fits', /deredden, aperture=150.0, $
      /meanprofile, /tracespec, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits

    ispec, 'ngc_0520_drift_100.fits', /deredden, aperture=155.0, $
      /meanprofile, /tracespec, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits

    ispec, 'ngc_0523_drift_030.fits', /deredden, refrow=75L, aperture=100.0, $
      /meanprofile, tracespec=0, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits

; the spectrum is centered on NGC0750, not NGC0751; spatial overlap 
    ispec, 'ngc_0750_drift_030.fits', /deredden, aperture=120.0, $
      /meanprofile, /tracespec, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits

; ###########################################################################
; NGC0943 is at row ~120 and NGC0942 is at row ~80; significant
; spatial overlap

; NGC0942
    ispec, 'ngc_0942-3_drift_033.fits', /deredden, refrow=63L, aperture=60.0, $
      /meanprofile, tracespec=0, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, telluric=telluric, $
      outname='ngc_0942_drift_033', galaxy='NGC0942'

; NGC0943
    ispec, 'ngc_0942-3_drift_033.fits', /deredden, refrow=135L, aperture=60.0, $
      /meanprofile, tracespec=0, traceorder=1, s1, outname='ngc_0943_drift_033', galaxy='NGC0943', $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, telluric=telluric

; ARP309 = NGC0942-3
    ispec, 'ngc_0942-3_drift_033.fits', /deredden, refrow=100L, aperture=130.0, $
      /meanprofile, tracespec=0, traceorder=1, s1, outname='arp_309_drift_033', galaxy='ARP309', $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, telluric=telluric
; ###########################################################################

    ispec, 'ngc_1614_drift_025.fits', /deredden, aperture=80.0, $
      /meanprofile, /tracespec, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, telluric=telluric

    ispec, 'ngc_2623_drift_040.fits', /deredden, aperture=90.0, $
      /meanprofile, /tracespec, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, telluric=telluric

    ispec, 'ngc_2782_drift_060.fits', /deredden, aperture=100.0, $
      /meanprofile, /tracespec, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, telluric=telluric

    ispec, 'ngc_3303_drift_030.fits', /deredden, aperture=70.0, $
      /meanprofile, /tracespec, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, telluric=telluric
    
    ispec, 'ngc_3448_drift_060.fits', /deredden, aperture=160.0, $
      /meanprofile, /tracespec, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, telluric=telluric
    
    ispec, 'ngc_3656_drift_040.fits', /deredden, aperture=110.0, $
      /meanprofile, /tracespec, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, telluric=telluric
    
    ispec, 'ngc_3921_drift_020.fits', /deredden, aperture=100.0, $
      /meanprofile, /tracespec, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, telluric=telluric
 
; NGC4038-9 = ARP244 = Antennae; extract together; unable to extract
; individual spectra
    ispec, 'ngc_4038-9_drift_165.fits', /deredden, aperture=200.0, $
      /meanprofile, /tracespec, traceorder=1, s1, outname='arp_244_drift_165', galaxy='ARP244', $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, telluric=telluric

    ispec, 'ngc_4194_drift_030.fits', /deredden, aperture=125.0, $
      /meanprofile, /tracespec, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, telluric=telluric

; ###########################################################################
; two drift-scans were taken, one at 90 degrees with significant
; spatial overlap, and one at 180 degrees, which allows separation of
; the two individual galaxies; NGC4676 = ARP242

; POSANGLE=90; NGC4676A
    ispec, 'ngc_4676_drift_150.fits', /deredden, refrow=130, loaperture=8.0, upaperture=32.0, $
      /meanprofile, /tracespec, traceorder=1, sbox=10, s1, outname='ngc_4676a_drift_150', galaxy='NGC4676A', $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, telluric=telluric

; POSANGLE=90; NGC4676B    
    ispec, 'ngc_4676_drift_150.fits', /deredden, refrow=100, loaperture=30.0, upaperture=5.0, $
      /meanprofile, /tracespec, traceorder=1, sbox=10, s1, outname='ngc_4676b_drift_150', galaxy='NGC4676B', $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, telluric=telluric

; POSANGLE=90; NGC4676 (together)
    ispec, 'ngc_4676_drift_150.fits', /deredden, loaperture=43.0, upaperture=32.0, $
      /meanprofile, /tracespec, traceorder=1, sbox=10, s1, outname='ngc_4676_drift_150', galaxy='NGC4676', $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, telluric=telluric

; ##########    

; POSANGLE=180; NGC4676A
    ispec, 'ngc_4676_drift_155.fits', /deredden, refrow=150, loaperture=15.0, upaperture=30.0, $
      /meanprofile, /tracespec, traceorder=1, sbox=10, s2, outname='ngc_4676a_drift_155', galaxy='NGC4676A', $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, telluric=telluric

; POSANGLE=180; NGC4676B    
    ispec, 'ngc_4676_drift_155.fits', /deredden, refrow=110, loaperture=30.0, upaperture=20.0, $
      /meanprofile, /tracespec, traceorder=1, sbox=10, s2, outname='ngc_4676b_drift_155', galaxy='NGC4676B', $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, telluric=telluric

; POSANGLE=180; NGC4676 (together)
    ispec, 'ngc_4676_drift_155.fits', /deredden, loaperture=30.0, upaperture=65.0, $
      /meanprofile, /tracespec, traceorder=1, s2, outname='ngc_4676_drift_155', galaxy='NGC4676', $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, telluric=telluric

; ###########################################################################

; ###########################################################################
; NGC5278-9=ARP239; interacting pair; extract separately and together;
; significant spatial overlap 
    
; NGC5278
    ispec, 'ngc_5278-9_drift_042.fits', /deredden, refrow=155, loaperture=25.0, upaperture=30.0, $
      /meanprofile, /tracespec, traceorder=1, s1, outname='ngc_5278_drift_042', galaxy='NGC5278', $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, telluric=telluric

; NGC5279
    ispec, 'ngc_5278-9_drift_042.fits', /deredden, refrow=110, loaperture=20.0, upaperture=12.0, $
      /meanprofile, /tracespec, traceorder=1, sbox=10L, s1, outname='ngc_5279_drift_042', galaxy='NGC5279', $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, telluric=telluric

; ARP239 (together)
    ispec, 'ngc_5278-9_drift_042.fits', /deredden, loaperture=60.0, upaperture=30.0, $
      /meanprofile, /tracespec, traceorder=1, s1, outname='arp_239_drift_042', galaxy='ARP239', $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, telluric=telluric
; ###########################################################################

; ###########################################################################
; NGC7284-5=ARP093; interacting pair; significant spatial overlap

; NGC7284
    ispec, 'ngc_7284-5_drift_070.fits', /deredden, refrow=130, loaperture=7.0, upaperture=43.0, $
      /meanprofile, /tracespec, traceorder=1, sbox=10, s1, outname='ngc_7284_drift_070', galaxy='NGC7284', $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, telluric=telluric

; NGC7285    
    ispec, 'ngc_7284-5_drift_070.fits', /deredden, refrow=110, loaperture=44.0, upaperture=6.0, $
      /meanprofile, /tracespec, traceorder=1, s1, outname='ngc_7285_drift_070', galaxy='NGC7285', $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, telluric=telluric

; ARP093 (together)      
    ispec, 'ngc_7284-5_drift_070.fits', /deredden, loaperture=45.0, upaperture=55.0, $
      /meanprofile, /tracespec, traceorder=1, s1, outname='arp_093_drift_070', galaxy='ARP093', $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, telluric=telluric
; ###########################################################################

; ARP223    
    ispec, 'ngc_7585_drift_060.fits', /deredden, aperture=130.0, $
      /meanprofile, /tracespec, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, telluric=telluric

; ARP222    
    ispec, 'ngc_7727_drift_100.fits', /deredden, aperture=150.0, $
      /meanprofile, /tracespec, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, telluric=telluric

; ARP195; galaxy triple; extract together
    ispec, 'ugc_04653_drift_084.fits', /deredden, loaperture=50.0, upaperture=34.0, $
      /meanprofile, /tracespec, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, telluric=telluric

; ###########################################################################
; the SE component is UGC09425_NED02; the NW component is
; UGC09425_NED01; the composite system is called UGC09425=ARP241;
; spatial overlap

    ispec, 'ugc_09425_drift_025.fits', /deredden, refrow=140, loaperture=8.0, upaperture=27.0, $
      /meanprofile, /tracespec, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, telluric=telluric, $
      outname='ugc_09425nw_drift_025', galaxy='UGC09425NW'
    
    ispec, 'ugc_09425_drift_025.fits', /deredden, refrow=115, loaperture=27.0, upaperture=8.0, $
      /meanprofile, /tracespec, traceorder=1, sbox=10, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, telluric=telluric, $
      outname='ugc_09425se_drift_025', galaxy='UGC09425SE'
    
    ispec, 'ugc_09425_drift_025.fits', /deredden, loaperture=45.0, upaperture=25.0, $
      /meanprofile, /tracespec, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, telluric=telluric, $
      outname='ugc_09425_drift_025', galaxy='UGC09425'
; ###########################################################################
    
; ###########################################################################
; UGC11175=ARP081=NGC6621-2; NGC6621 is NW of NGC6622; significant
; spatial overlap

; NGC6621; *NB* include the tidal debris feature in this spectrum  
    ispec, 'ugc_11175_drift_060.fits', /deredden, refrow=120, loaperture=27.0, upaperture=58.0, $
      /meanprofile, /tracespec, traceorder=1, s1, outname='ngc_6621_drift_060', galaxy='NGC6621', $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, telluric=telluric

; NGC6622
    ispec, 'ugc_11175_drift_060.fits', /deredden, refrow=70, loaperture=27.0, upaperture=8.0, $
      /meanprofile, /tracespec, traceorder=1, s1, outname='ngc_6622_drift_060', galaxy='NGC6622', $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, telluric=telluric

; UGC11175 (together)
    ispec, 'ugc_11175_drift_060.fits', /deredden, refrow=120, loaperture=62.0, upaperture=58.0, $
      /meanprofile, /tracespec, traceorder=1, s1, outname='ugc_11175_drift_060', galaxy='UGC11175', $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, telluric=telluric
; ###########################################################################

stop    
    
return
end
    
