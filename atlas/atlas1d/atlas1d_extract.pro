pro atlas1d_extract, wfits=wfits, noplot=noplot
; J. Moustakas, 2002-2003, U of A
; (also see IUE_EXTRACT)

    datapath = atlas_path(/atlas2d)
    outpath = atlas_path(/atlas1d)
    repeatpath = atlas_path(/atlas1d)+'repeaters/'

    t0 = systime(1)

; ###########################################################################
; ARP256 interacting pair; 180 degree position angle; MCG-02-01-051 is
; south of MCG-02-01-052 (higher row number)
    
; ARP256NED01; MCG-02-01-051
    ispec, 'arp_256_drift_040.fits', /deredden, refrow=82, aperture=40, $
      /meanprofile, tracespec=1, traceorder=1, fluxrange=[-10,200]*1E-17, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='arp_256s_drift_040', galaxy='ARP256S'

; very small amount of stellar contamination that I am leaving out of
; the aperture; ARP256NED02; MCG-02-01-052
    ispec, 'arp_256_drift_040.fits', /deredden, refrow=51, loaperture=22, upaperture=38, $
      /meanprofile, tracespec=1, traceorder=1, sbox=10, fluxrange=[-10,200]*1E-17, s1, $
;     /starfit, star_refrow=28, debug=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='arp_256n_drift_040', galaxy='ARP256N'

; ARP256 (together)
    ispec, 'arp_256_drift_040.fits', /deredden, refrow=82, loaperture=80.0, upaperture=20, $
      /meanprofile, tracespec=1, traceorder=1, fluxrange=[-10,200]*1E-17, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='arp_256_drift_040', galaxy='ARP256'
; ###########################################################################

    ispec, 'cgcg_049-057_drift_030.fits', /deredden, aperture=35, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; ###########################################################################
; interacting pair (the pair is also called CGCG 239-011); 90 degree
; position angle; MCG+08-18-012 is West (higher row) number relative
; to MCG+08-18-013

; CGCG239-011NED01; WEST; MCG+08-18-012
    ispec, 'cgcg_239-011_drift_030.fits', /deredden, refrow=62, aperture=50, $
      /meanprofile, tracespec=1, traceorder=1, fluxrange=[-0.5,4]*1E-16, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='cgcg_239-011w_drift_030', galaxy='CGCG239-011W'
    
; CGCG239-011NED02; EAST; MCG+08-18-013
    ispec, 'cgcg_239-011_drift_030.fits', /deredden, refrow=22, loaperture=15, upaperture=25, $
      /meanprofile, tracespec=1, traceorder=1, sbox=10, fluxrange=[-0.1,2]*1E-16, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='cgcg_239-011e_drift_030', galaxy='CGCG239-011E'
; ###########################################################################
    
    ispec, 'cgcg_436-030_drift_040.fits', /deredden, refrow=69, aperture=35, $
      /meanprofile, tracespec=1, traceorder=1, sbox=10, fluxrange=[-0.2,1.5]*1E-15, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; dusty!
    ispec, 'cgcg_453-062_drift_030.fits', /deredden, aperture=50, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; ###########################################################################
; interacting pair; MCG-03-12-002; 180 degree position angle; extract
; both components and the N and S components separately; small amount
; of stellar contamination in the summed and the N components; the N
; component is at smaller row number; the S component is of lower S/N;
; setting MEANPROFILE is less instructive here

; together; foreground star subtracted
    ispec, 'eso_550-ig_025_drift_040.fits', /deredden, refrow=60, loaperture=20.0, upaperture=32.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      /starfit, star_refrow=46, trace_searchbox=4L, frac_rows=0.04, debug=debug, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      galaxy='ESO550-IG025'
    
; N (ESO550-IG025 NED01); foreground star subtracted
    ispec, 'eso_550-ig_025_drift_040.fits', /deredden, loaperture=20.0, upaperture=12.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      /starfit, star_refrow=46, trace_searchbox=4L, frac_rows=0.04, debug=debug, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='eso_550-ig_025n_drift_040', galaxy='ESO550-IG025N'

; S (ESO550-IG025 NED02)
    ispec, 'eso_550-ig_025_drift_040.fits', /deredden, refrow=70, loaperture=7, upaperture=13, $
      /meanprofile, tracespec=1, traceorder=1, sbox=5, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='eso_550-ig_025s_drift_040', galaxy='ESO550-IG025S'
; ###########################################################################

; MCG-03-57-017; possible CR on top of OIII 4959
    ispec, 'eso_602-g_025_drift_030.fits', /deredden, aperture=90, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric
    
; MCG+01-31-030; CGCG041-077
    ispec, 'haro_06_drift_020.fits', /deredden, aperture=30, $
      /meanprofile, tracespec=1, traceorder=2, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; GTO starburst
    ispec, 'hs_0822+3542_drift_020.fits', /deredden, refrow=133, aperture=20, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; ARP230; the object on the shoulder is a residual cosmic ray; the
; trace fails, so do not trace; Turner
    ispec, 'ic_0051_drift_089.fits', /deredden, aperture=75.0, $
      /meanprofile, tracespec=0, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits

; ARP228 = UGC01267; Turner
    ispec, 'ic_0162_drift_040.fits', /deredden, aperture=100.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits

    ispec, 'ic_0691_drift_040.fits', /deredden, aperture=40, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ic_0749_drift_090.fits', /deredden, aperture=140.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; dusty!    
    ispec, 'ic_0750_drift_090.fits', /deredden, aperture=100, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric
    
    ispec, 'ic_0860_drift_030.fits', /deredden, aperture=40, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; ###########################################################################
; repeat observations; observed 5/99
    ispec, 'ic_0883_drift_030.fits', /deredden, aperture=50.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; Turner merger
    ispec, 'ic_0883_drift_1_020.fits', /deredden, aperture=90.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits

    ispec, 'ic_0883_drift_2_020.fits', /deredden, aperture=90.0, $
      /meanprofile, tracespec=1, traceorder=1, s2, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits
; ###########################################################################

    ispec, 'ic_1076_drift_060.fits', /deredden, aperture=40, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ic_1586_drift_020.fits', /deredden, aperture=30, $
      /meanprofile, tracespec=1, traceorder=2, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; ###########################################################################
; interacting pair; IC1623A/IC1623B (the system is called IC1623);
; extract separately and together; 90 degree position angle; difficult
; extraction  

; IC1623A
    ispec, 'ic_1623_drift_030.fits', /deredden, refrow=59, aperture=30, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='ic_1623a_drift_030', galaxy='IC1623A'

; IC1623B
    ispec, 'ic_1623_drift_030.fits', /deredden, refrow=37, aperture=20, $
      /meanprofile, tracespec=1, traceorder=1, sbox=7, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='ic_1623b_drift_030', galaxy='IC1623B'

; IC1623 = ARP236 (together)
    ispec, 'ic_1623_drift_030.fits', /deredden, refrow=59, loaperture=35, upaperture=15.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='ic_1623_drift_030', galaxy='IC1623'
; ###########################################################################

; foreground stellar contamination    
    ispec, 'ic_1727_drift_240.fits', /deredden, refrow=65, aperture=180, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; ###########################################################################
; interacting pair; also called UGC06436; IC2810B (South-East of
; IC2810) has very low S/N compared to IC2810; exclude IC2810B from
; the atlas (jm05jul27uofa)

    ispec, 'ic_2810_drift_060.fits', /deredden, refrow=32, aperture=30, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=0, /noskymask, telluric=telluric, $
      outname='ic_2810b_drift_060', galaxy='IC2810B'
    
    ispec, 'ic_2810_drift_060.fits', /deredden, refrow=75, loaperture=15.0, upaperture=20.0, $
      /meanprofile, tracespec=1, traceorder=1, sbox=10, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='ic_2810_drift_060', galaxy='IC2810'
; ###########################################################################

; ARP220    
    ispec, 'ic_4553_drift_036.fits', /deredden, aperture=80.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ic_5179_drift_060.fits', /deredden, aperture=130, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ic_5298_drift_030.fits', /deredden, aperture=40.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; interacting pair; extract together
    ispec, 'iii_zw_035_drift_020.fits', /deredden, refrow=57, aperture=35, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; ###########################################################################
; NW component    
    ispec, 'ii_zw_096_drift_040.fits', /deredden, loaperture=10.0, upaperture=20.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='ii_zw_096nw_drift_040', galaxy='IIZW096NW'

; SE component    
    ispec, 'ii_zw_096_drift_040.fits', /deredden, refrow=46, loaperture=12.0, upaperture=6.0, $
      /meanprofile, tracespec=1, traceorder=1, sbox=5, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='ii_zw_096se_drift_040', galaxy='IIZW096SE'

; together    
    ispec, 'ii_zw_096_drift_040.fits', /deredden, loaperture=30.0, upaperture=20.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      galaxy='IIZW096'
; ###########################################################################

; interacting pair (2MASXJ03384633+1532549, 2MASXJ03384713+1532537);
; extract together; unable to get two good spectra
    ispec, 'iras_03359+1523_drift_015.fits', /deredden, refrow=68, aperture=35, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; AGN
    ispec, 'iras_05189-2524_drift_030.fits', /deredden, aperture=30, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; 2MASXIJ0900253+390354; interacting pair
    ispec, 'iras_08572+3915_drift_020.fits', /deredden, refrow=58, aperture=30, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; IRAS15250+3609 has been excluded from the atlas (02apr) because
; H-alpha is redshifted beyond the spectral coverage
    ispec, 'iras_15250+3609_drift_030.fits', /deredden, aperture=15.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=0, /noskymask, telluric=telluric

; interacting pair
    ispec, 'iras_17132+5313_drift_020.fits', /deredden, refrow=57, aperture=30, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; very dusty!
    ispec, 'iras_17208-0014_drift_015.fits', /deredden, aperture=25.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; 2MASXIJ2339012+362108
    ispec, 'iras_23365+3604_drift_020.fits', /deredden, aperture=25.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; interacting pair; IRAS12596-1529
    ispec, 'mcg-02-33-098_drift_050.fits', /deredden, aperture=90.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric
    
; (also called IRAS 01077-1707 and 2MASXIJ0110090-165109)
    ispec, 'mcg-03-04-014_drift_020.fits', /deredden, aperture=35.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric    

; IUE galaxy; interacting pair CGCG288-028=KUG0912+599; MRK019 =
; CGCG288-028NED01 is NE of CGCG288-028NED02; the drift-scan aperture
; only includes MRK0019
    ispec, 'mrk_0019_drift_015.fits', /deredden, aperture=40, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'mrk_0026_drift_020.fits', /deredden, aperture=25, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'mrk_0053_drift_020.fits', /deredden, aperture=30.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; ###########################################################################
; interacting pair; MRK0055 = PGC044300 is W and
; 2MASSXJ12572600+2724188 = PGC1807147 is E; the E component is low
; S/N, so exclude from the atlas
    ispec, 'mrk_0055_drift_020.fits', /deredden, loaperture=5.0, upaperture=15.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      galaxy='MRK0055'

    ispec, 'mrk_0055_drift_020.fits', /deredden, refrow=55, loaperture=10.0, upaperture=5.0, $
      /meanprofile, tracespec=1, traceorder=1, sbox=5, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=0, /noskymask, telluric=telluric, $
      outname='mrk_0055e_drift_020', galaxy='MRK0055E'
; ###########################################################################

    ispec, 'mrk_0066_drift_020.fits', /deredden, refrow=67, aperture=30, $
      /meanprofile, tracespec=1, traceorder=1, sbox=5, fluxrange=[-2,10]*1E-16, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'mrk_0247_drift_014.fits', /deredden, aperture=30, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'mrk_0296_drift_040.fits', /deredden, aperture=30, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'mrk_0312_drift_020.fits', /deredden, loaperture=17.0, upaperture=13.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'mrk_0315_drift_020.fits', /deredden, aperture=35, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; ###########################################################################
; GTO starburst; the first scan has higher S/N 
    ispec, 'mrk_0331_drift_1_040.fits', /deredden, aperture=40.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'mrk_0331_drift_2_040.fits', /deredden, aperture=40.0, $
      /meanprofile, tracespec=1, traceorder=1, s2, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric
; ###########################################################################

    ispec, 'mrk_0360_drift_030.fits', /deredden, aperture=40.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'mrk_0475_drift_020.fits', /deredden, aperture=30.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric
    
; IUE galaxy (used to be UGC10565)
    ispec, 'mrk_0499_drift_020.fits', /deredden, aperture=30, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'mrk_0600_drift_030.fits', /deredden, loaperture=13.0, upaperture=15.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'mrk_0685_drift_020.fits', /deredden, aperture=40.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; interacting pair with significant spatial overlap; 180-degree slit
; position angle, although the position of the star is confusing;
; assuming 180, North is toward lower row numbers; the N component has
; very low signal; just extract the sum, which saves me from having to
; identify (possibly incorrectly) the N and S components
    ispec, 'mrk_0848_drift_030.fits', /deredden, refrow=63, loaperture=35.0, upaperture=15.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'mrk_0863_drift_020.fits', /deredden, refrow=65, aperture=50, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; GTO starburst
    ispec, 'mrk_0930_drift_015.fits', /deredden, loaperture=15.0, upaperture=10.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'mrk_0960_drift_040.fits', /deredden, loaperture=20.0, upaperture=20.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; GTO starburst
    ispec, 'mrk_1450_drift_015.fits', /deredden, aperture=20.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'mrk_1460_drift_020.fits', /deredden, aperture=20.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'mrk_1490_drift_020.fits', /deredden, aperture=40, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric
    
; subtracted supernova SN1955C
    ispec, 'ngc_0023_drift_060.fits', /deredden, refrow=63, aperture=70, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      /starfit, star_refrow=49, frac_rows=0.06, debug=debug, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_0034_drift_040.fits', /deredden, aperture=45, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_0095_drift_060.fits', /deredden, aperture=100, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_0151_drift_150.fits', /deredden, aperture=120, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; foreground star excluded
    ispec, 'ngc_0157_drift_180.fits', /deredden, refrow=49, loaperture=50.0, upaperture=80.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_0178_drift_090.fits', /deredden, aperture=60, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; dusty!    
    ispec, 'ngc_0232_drift_040.fits', /deredden, aperture=60.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_0244_drift_060.fits', /deredden, aperture=70, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric 

    ispec, 'ngc_0245_drift_090.fits', /deredden, aperture=80, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_0278_drift_090.fits', /deredden, aperture=110, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_0309_drift_120.fits', /deredden, refrow=63, aperture=150, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; SINGS galaxy    
    ispec, 'ngc_0337_drift_120.fits', /deredden, aperture=150, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_0474_drift_060.fits', /deredden, aperture=150.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits

    ispec, 'ngc_0520_drift_100.fits', /deredden, aperture=155.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits

    ispec, 'ngc_0523_drift_030.fits', /deredden, refrow=75L, aperture=100.0, $
      /meanprofile, tracespec=0, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits

    ispec, 'ngc_0615_drift_150.fits', /deredden, aperture=90.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; dusty!
    ispec, 'ngc_0660_drift_120.fits', /deredden, refrow=62, loaperture=60.0, upaperture=70.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_0694_drift_030.fits', /deredden, aperture=40.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_0695_drift_030.fits', /deredden, refrow=59, aperture=45, $
      /meanprofile, tracespec=1, traceorder=1, sbox=10, fluxrange=[-0.1,2]*1E-15, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; the spectrum is centered on NGC0750, not NGC0751; spatial overlap 
    ispec, 'ngc_0750_drift_030.fits', /deredden, aperture=80.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits

    ispec, 'ngc_0784_drift_210.fits', /deredden, refrow=62, aperture=80, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_0877_drift_090.fits', /deredden, aperture=120, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_0922_drift_120.fits', /deredden, aperture=100.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; ###########################################################################
; Turner; NGC0943 is at row ~120 and NGC0942 is at row ~80;
; significant spatial overlap

; NGC0942
    ispec, 'ngc_0942-3_drift_033.fits', /deredden, refrow=82, loaperture=40.0, upaperture=15.0, $
      /meanprofile, tracespec=1, traceorder=1, sbox=5, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='ngc_0942_drift_033', galaxy='NGC0942'

; NGC0943
    ispec, 'ngc_0942-3_drift_033.fits', /deredden, refrow=120, loaperture=15.0, upaperture=40.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='ngc_0943_drift_033', galaxy='NGC0943'

; ARP309 = NGC0942-3
    ispec, 'ngc_0942-3_drift_033.fits', /deredden, loaperture=70.0, upaperture=40.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='arp_309_drift_033', galaxy='ARP309'
; ###########################################################################

    ispec, 'ngc_0958_drift_120.fits', /deredden, aperture=60.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_0959_drift_060.fits', /deredden, aperture=130.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_0972_drift_090.fits', /deredden, aperture=110.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_0976_drift_060.fits', /deredden, aperture=80, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric
    
; ###########################################################################
; repeat observations; the second scan has much higher S/N);
; foreground star excluded

    ispec, 'ngc_1003_drift_1_160.fits', /deredden, refrow=63, loaperture=20.0, upaperture=50.0, $
      /meanprofile, tracespec=1, traceorder=1, sbox=10, fluxrange=[-30,170]*1E-17, s1, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_1003_drift_2_160.fits', /deredden, refrow=63, loaperture=20.0, upaperture=50.0, $
      /meanprofile, tracespec=1, traceorder=1, sbox=10, fluxrange=[-30,500]*1E-17, s2, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric
; ###########################################################################
    
    ispec, 'ngc_1036_drift_060.fits', /deredden, aperture=60.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; foreground star subtracted
    ispec, 'ngc_1058_drift_075.fits', /deredden, refrow=58, aperture=170, $
      /meanprofile, tracespec=1, traceorder=1, sbox=10, s1, $
      /starfit, star_refrow=72, debug=debug, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; IUE galaxy
    ispec, 'ngc_1068_drift_180.fits', /deredden, refrow=57, loaperture=75.0, upaperture=80.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_1084_drift_180.fits', /deredden, aperture=140, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_1087_drift_180.fits', /deredden, aperture=140, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; NGC1097 --> IUE_EXTRACT

    ispec, 'ngc_1134_drift_060.fits', /deredden, refrow=63, aperture=75, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; GTO starburst    
    ispec, 'ngc_1140_drift_120.fits', /deredden, aperture=60.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; ###########################################################################
; interacting pair; 90-degree slit position angle; blended; also
; called NGC1141-2; NGC1143 is West (higher row number) of NGC1144 

; NGC1143    
    ispec, 'ngc_1143-4_drift_060.fits', /deredden, refrow=65, loaperture=15.0, upaperture=30.0, $
      /meanprofile, tracespec=1, traceorder=1, sbox=5, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='ngc_1143_drift_060', galaxy='NGC1143'

; NGC1144
    ispec, 'ngc_1143-4_drift_060.fits', /deredden, refrow=42, aperture=50.0, $
      /meanprofile, tracespec=1, traceorder=1, s2, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='ngc_1144_drift_060', galaxy='NGC1144'

; ARP 118 (together)
    ispec, 'ngc_1143-4_drift_060.fits', /deredden, refrow=42, loaperture=25.0, upaperture=70.0, $
      /meanprofile, tracespec=1, traceorder=1, s3, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='arp_118_drift_060', galaxy='ARP118'
; ###########################################################################
    
; foreground stellar contamination? no, I think the points seen in the
; DSS image are actually HII regions, determined by looking at an
; H-alpha image
    ispec, 'ngc_1156_drift_080.fits', /deredden, refrow=62, aperture=160.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; SINGS
    ispec, 'ngc_1266_drift_055.fits', /deredden, refrow=63, aperture=90, $
      /meanprofile, tracespec=1, traceorder=1, fluxrange=[-10,140]*1E-17, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; AGN!
    ispec, 'ngc_1275_drift_040.fits', /deredden, aperture=75, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric
    
    ispec, 'ngc_1345_drift_090.fits', /deredden, aperture=80.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_1359_drift_090.fits', /deredden, refrow=67, loaperture=40.0, upaperture=70.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; SINGS
    ispec, 'ngc_1377_drift_055.fits', /deredden, aperture=110, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric
    
    ispec, 'ngc_1385_drift_180.fits', /deredden, refrow=61, aperture=110, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_1421_drift_240.fits', /deredden, refrow=50, loaperture=38.0, upaperture=27.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric
    
; foreground stellar contamination; foreground star subtracted 
    ispec, 'ngc_1560_drift_600.fits', /deredden, refrow=65, aperture=90, $
      /meanprofile, tracespec=1, traceorder=1, fluxrange=[-100,1300]*1E-17, s1, $
      /starfit, star_refrow=67, debug=debug, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; IUE galaxy; foreground stellar contamination
    ispec, 'ngc_1569_drift_060.fits', /deredden, refrow=57, aperture=130, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; ###########################################################################
; repeat observations
    ispec, 'ngc_1614_drift_025.fits', /deredden, aperture=80.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; IUE galaxy; the second scan has higher S/N; GTO starburst
    ispec, 'ngc_1614_drift_1_060.fits', /deredden, aperture=80, $
      /meanprofile, tracespec=1, traceorder=1, s2, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_1614_drift_2_060.fits', /deredden, aperture=80, $
      /meanprofile, tracespec=1, traceorder=1, s3, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric
; ###########################################################################
      
; IUE galaxy
    ispec, 'ngc_1800_drift_060.fits', /deredden, refrow=52, aperture=100, $
      /meanprofile, tracespec=1, traceorder=1, fluxrange=[-0.7,3]*1E-15, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_2139_drift_120.fits', /deredden, aperture=120, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; dusty!  foreground stellar contamination; GTO starburst
    ispec, 'ngc_2146_drift_180.fits', /deredden, aperture=200, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric
    
    ispec, 'ngc_2337_drift_050.fits', /deredden, refrow=65, aperture=90, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; ###########################################################################
    ispec, 'ngc_2363a_drift_090.fits', /deredden, refrow=110, aperture=30.0, $
      /meanprofile, tracespec=1, traceorder=1, fluxrange=[-0.1,0.9]*1E-15, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='ngc_2363_drift_090', galaxy='NGC2363'

    ispec, 'ngc_2363a_drift_090.fits', /deredden, loaperture=40.0, upaperture=20.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      galaxy='NGC2363A'
; ###########################################################################

 ; dusty!    
    ispec, 'ngc_2388_drift_030.fits', /deredden, aperture=60.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric
    
; IUE galaxy
    ispec, 'ngc_2415_drift_050.fits', /deredden, aperture=55, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; foreground stellar contamination
    ispec, 'ngc_2500_drift_120.fits', /deredden, refrow=68, loaperture=65.0, upaperture=85.0, $
      /meanprofile, tracespec=1, traceorder=1, fluxrange=[-0.5,3]*1E-15, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; IUE galaxy; GTO starburst
    ispec, 'ngc_2537_drift_060.fits', /deredden, refrow=57, aperture=100, $
      /meanprofile, tracespec=1, traceorder=1, sbox=10, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; foreground stellar contamination    
    ispec, 'ngc_2541_drift_240.fits', /deredden, refrow=63, aperture=170, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_2552_drift_120.fits', /deredden, aperture=170.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_2599_drift_040.fits', /deredden, aperture=100.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_2623_drift_040.fits', /deredden, aperture=90.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; ###########################################################################
; repeat observations; GTO starburst

; Turner    
    ispec, 'ngc_2782_drift_060.fits', /deredden, aperture=110.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_2782_drift_100.fits', /deredden, aperture=110, $
      /meanprofile, tracespec=1, traceorder=1, s2, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric
; ###########################################################################

    ispec, 'ngc_2893_drift_040.fits', /deredden, aperture=60.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; GTO starburst
    ispec, 'ngc_2903_drift_360.fits', /deredden, aperture=200, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; IUE galaxy
    ispec, 'ngc_3049_drift_120.fits', /deredden, aperture=100.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; GTO starburst
    ispec, 'ngc_3077_drift_200.fits', /deredden, aperture=200, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; ###########################################################################
; the second scan is better positioned spatially and has higher S/N;
; stellar contamination; exclude the stars near row 90 from the
; aperture; difficult extraction; GTO starburst
    ispec, 'ngc_3079_drift_1_360.fits', /deredden, refrow=80, loaperture=40.0, upaperture=50.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_3079_drift_2_360.fits', /deredden, refrow=57, loaperture=50.0, upaperture=50.0, $
      /meanprofile, tracespec=1, traceorder=1, s2, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric
; ###########################################################################

; foreground star subtracted
    ispec, 'ngc_3104_drift_105.fits', /deredden, refrow=65, loaperture=80.0, upaperture=60.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, inputspec=inputspec, $
      /starfit, star_refrow=68, debug=debug, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_3198_drift_480.fits', /deredden, aperture=220.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric
    
; foreground star excluded; damn that's bright!; attempted to subtract
; (jm05jul20uofa, but the star is just too bright)
    ispec, 'ngc_3239_drift_090.fits', /deredden, refrow=60, loaperture=50, upaperture=40.0, $
      /meanprofile, tracespec=1, traceorder=1, sbox=5, fluxrange=[-1,6]*1E-15, s1, $      
;     /starfit, star_refrow=90, debug=debug, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; SINGS; foreground star excluded
    ispec, 'ngc_3265_drift_055.fits', /deredden, aperture=45, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_3274_drift_040.fits', /deredden, aperture=120.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; Turner    
    ispec, 'ngc_3303_drift_030.fits', /deredden, aperture=70.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_3310_drift_065.fits', /deredden, aperture=90.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; stitched together; three foreground stars subtracted; the star at
; row 58 is not subtracted perfectly
    ispec, 'ngc_3344_drift_360.fits', /deredden, refrow=120, aperture=375.0, $
      /meanprofile, tracespec=1, traceorder=1, sbox=10, fluxrange=[-40,700]*1E-17, s2, $
      /starfit, star_refrow=[58,82,105], frac_rows=0.04, debug=debug, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_3351_drift_360.fits', /deredden, aperture=200, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_3353_drift_060.fits', /deredden, aperture=70, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_3365_drift_180.fits', /deredden, aperture=110, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; GTO starburst    
    ispec, 'ngc_3367_drift_100.fits', /deredden, aperture=140, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; ###########################################################################
; NGC3395/6 = ARP 270; both NGC3395 and NGC3596 are IUE galaxies;
; 90-degree slit position angle; NGC3395 is West of NGC3396

; NGC3395 = ARP 270 W    
    ispec, 'ngc_3395-6_drift_080.fits', /deredden, refrow=89, loaperture=30.0, upaperture=50.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='ngc_3395_drift_080', galaxy='NGC3395'

; NGC3396 = ARP 270 E
    ispec, 'ngc_3395-6_drift_080.fits', /deredden, refrow=45, loaperture=52.0, upaperture=38.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='ngc_3396_drift_080', galaxy='NGC3396'

; ARP 270 (together)
    ispec, 'ngc_3395-6_drift_080.fits', /deredden, refrow=89, loaperture=120.0, upaperture=50.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='arp_270_drift_080', galaxy='ARP270'
; ###########################################################################

; IUE galaxy; foreground stars excluded
    ispec, 'ngc_3432_drift_300.fits', /deredden, loaperture=35.0, upaperture=20.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_3442_drift_030.fits', /deredden, aperture=40, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; ###########################################################################
; IUE galaxy

; Turner
    ispec, 'ngc_3448_drift_060.fits', /deredden, aperture=160.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_3448_drift_090.fits', /deredden, aperture=140, $
      /meanprofile, tracespec=1, traceorder=1, s2, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric
; ###########################################################################

; foreground star subtracted
    ispec, 'ngc_3504_drift_120.fits', /deredden, aperture=90.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      /starfit, star_refrow=30, frac_rows=0.05, debug=debug, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_3510_drift_120.fits', /deredden, aperture=70.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_3521_drift_360.fits', /deredden, aperture=190.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_3600_drift_160.fits', /deredden, aperture=60.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; why is the S/N so high?  GTO starburst
    ispec, 'ngc_3628_drift_550.fits', /deredden, refrow=40.0, loaperture=60.0, upaperture=120.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; Turner    
    ispec, 'ngc_3656_drift_040.fits', /deredden, aperture=110.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric
    
; ###########################################################################
; IUE galaxy; interacting pair; 90 degree position angle; extract the
; sum (NGC3690), the North-East component (MCG+10-17-005) and the
; South-West component (MCG+10-17-003) separately 

; summed spectrum (NGC3690)    
    ispec, 'ngc_3690_drift_060.fits', /deredden, refrow=65, loaperture=55.0, upaperture=35.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      galaxy='NGC3690'

; North-East component (MCG+10-17-005, NGC3690NED02)
    ispec, 'ngc_3690_drift_060.fits', /deredden, refrow=52, loaperture=35.0, upaperture=10.0, $
      /meanprofile, tracespec=1, traceorder=1, sbox=10, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='ngc_3690ne_drift_060', galaxy='NGC3690NE'
    
; South-West component (MCG+10-17-003, NGC3690NED01)
    ispec, 'ngc_3690_drift_060.fits', /deredden, refrow=78, loaperture=10.0, upaperture=35.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='ngc_3690sw_drift_060', galaxy='NGC3690SW'
; ###########################################################################

    ispec, 'ngc_3718_drift_180.fits', /deredden, aperture=180, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric
    
    ispec, 'ngc_3726_drift_240.fits', /deredden, aperture=200, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric
    
; ###########################################################################
; repeat observations
    ispec, 'ngc_3729_drift_1_090.fits', /deredden, aperture=110.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric
    
    ispec, 'ngc_3729_drift_2_090.fits', /deredden, aperture=110.0, $
      /meanprofile, tracespec=1, traceorder=1, s2, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric    
; ###########################################################################

; IUE galaxy
    ispec, 'ngc_3738_drift_060.fits', /deredden, aperture=100.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; foreground stars excluded
    ispec, 'ngc_3741_drift_060.fits', /deredden, refrow=60, aperture=55.0, $
      /meanprofile, tracespec=1, traceorder=1, sbox=10, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; ###########################################################################
; interacting pair (NGC3769 with NGC3769A)

    ispec, 'ngc_3769_drift_120.fits', /deredden, aperture=90, $
      /meanprofile, tracespec=1, traceorder=1, sbox=10, fluxrange=[-10,200]*1E-17, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='ngc_3769_drift_120', galaxy='NGC3769'

    ispec, 'ngc_3769_drift_120.fits', /deredden, refrow=32, loaperture=20.0, upaperture=25.0, $
      /meanprofile, tracespec=1, traceorder=1, sbox=10, fluxrange=[-10,70]*1E-17, s2, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='ngc_3769a_drift_120', galaxy='NGC3769A'

; ARP 280 (together)
    ispec, 'ngc_3769_drift_120.fits', /deredden, loaperture=90.0, upaperture=45.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='arp_280_drift_120', galaxy='ARP280'
; ###########################################################################
    
; SINGS; foreground stellar contamination; unable to subtract the star
; well, so leave it in
    ispec, 'ngc_3773_drift_055.fits', /deredden, aperture=70.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
;     /starfit, star_refrow=72, frac_rows=0.05, trace_searchbox=3L, debug=debug, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; foreground star subtracted
    ispec, 'ngc_3782_drift_040.fits', /deredden, refrow=60, loaperture=45.0, upaperture=55.0, $
      /meanprofile, tracespec=1, traceorder=1, sbox=10, fluxrange=[-0.1,1.2]*1E-15, s1, $
      /starfit, star_refrow=82, frac_rows=0.06, debug=debug, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_3870_drift_040.fits', /deredden, aperture=50, $
      /meanprofile, tracespec=1, traceorder=1, sbox=10, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_3877_drift_300.fits', /deredden, aperture=140, $
      /meanprofile, tracespec=1, traceorder=1, sbox=10, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; foreground star subtracted
    ispec, 'ngc_3893_drift_180.fits', /deredden, refrow=63, aperture=170, $
      /meanprofile, tracespec=1, traceorder=1, s1, $ 
      /starfit, star_refrow=88, debug=debug, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; foreground star subtracted (M giant)
    ispec, 'ngc_3896_drift_040.fits', /deredden, aperture=80.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      /starfit, star_refrow=77, frac_rows=0.04, debug=debug, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_3906_drift_090.fits', /deredden, refrow=65, loaperture=50.0, upaperture=60.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_3913_drift_120.fits', /deredden, aperture=120.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; foreground star subtracted
    ispec, 'ngc_3917_drift_240.fits', /deredden, refrow=60, aperture=110.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      /starfit, star_refrow=85, frac_rows=0.06, debug=debug, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; Turner    
    ispec, 'ngc_3921_drift_020.fits', /deredden, aperture=90.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric
 
    ispec, 'ngc_3928_drift_060.fits', /deredden, aperture=90.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric
    
    ispec, 'ngc_3949_drift_090.fits', /deredden, aperture=130, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; foreground star subtracted
    ispec, 'ngc_3953_drift_300.fits', /deredden, aperture=170, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      /starfit, star_refrow=95, frac_rows=0.06, debug=debug, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_3972_drift_160.fits', /deredden, refrow=60, aperture=120, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; IUE galaxy
    ispec, 'ngc_3982_drift_090.fits', /deredden, aperture=110, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; ###########################################################################
; use NGC3985 to study aperture effects; comparable S/N
    ispec, 'ngc_3985_drift_050.fits', /deredden, aperture=80, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_3985_drift_040.fits', /deredden, aperture=80, $
      /meanprofile, tracespec=1, traceorder=1, s2, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric
; ###########################################################################

; ###########################################################################
; IUE galaxy; galaxy (interacting) pair

; NGC3991N    
    ispec, 'ngc_3991_drift_060.fits', /deredden, refrow=60, loaperture=26.0, upaperture=9.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='ngc_3991n_drift_060', galaxy='NGC3991N'

; NGC3991S
    ispec, 'ngc_3991_drift_060.fits', /deredden, refrow=70, loaperture=6.0, upaperture=30.0, $
      /meanprofile, tracespec=1, traceorder=1, sbox=5, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='ngc_3991s_drift_060', galaxy='NGC3991S'

; together    
    ispec, 'ngc_3991_drift_060.fits', /deredden, refrow=70, loaperture=25.0, upaperture=45.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='ngc_3991_drift_060', galaxy='NGC3991'
; ###########################################################################

; M109; do not trace; exclude from atlas; low S/N; also Rob's notes
; indicate "scan terminated by jump"
    ispec, 'ngc_3992_drift_600.fits', /deredden, refrow=65, aperture=180.0, $
      /meanprofile, tracespec=0, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=0, /noskymask, telluric=telluric

; IUE galaxy; ARP 313 W
    ispec, 'ngc_3994_drift_060.fits', /deredden, aperture=40.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; IUE galaxy; ARP 313 E
    ispec, 'ngc_3995_drift_120.fits', /deredden, refrow=60, aperture=90, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_3998_drift_120.fits', /deredden, aperture=130.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_4004_drift_090.fits', /deredden, refrow=62, loaperture=30.0, upaperture=20.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_4010_drift_240.fits', /deredden, refrow=70, loaperture=40.0, upaperture=60.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; ###########################################################################
; repeat observations
    ispec, 'ngc_4020_drift_1_120.fits', /deredden, aperture=60.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_4020_drift_2_120.fits', /deredden, aperture=60.0, $
      /meanprofile, tracespec=1, traceorder=1, s2, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric
; ###########################################################################

; Turner; NGC4038-9 = ARP244 = Antennae; extract together; unable to
; extract individual spectra
    ispec, 'ngc_4038-9_drift_165.fits', /deredden, aperture=200.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='arp_244_drift_165', galaxy='ARP244'

    ispec, 'ngc_4051_drift_240.fits', /deredden, aperture=200.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_4062_drift_170.fits', /deredden, aperture=160.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric
    
; two foreground stars subtracted; the star at row ~68 is a G-type
; star while the one at row ~78 is an M giant
    ispec, 'ngc_4068_drift_120.fits', /deredden, refrow=55, aperture=110.0, $
      /meanprofile, tracespec=1, traceorder=1, fluxrange=[-0.2,1.5]*1E-15, s1, $
      /starfit, star_refrow=[68,79], trace_searchbox=4L, frac_rows=0.04, debug=debug, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_4085_drift_040.fits', /deredden, aperture=160, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_4088_drift_300.fits', /deredden, aperture=140, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric
    
; ###########################################################################
; repeat observations; the first scan has very low S/N compared to the
; second scan (S/N = 7 versus 23)

    ispec, 'ngc_4096_drift_1_380.fits', /deredden, aperture=110, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_4096_drift_2_380.fits', /deredden, aperture=110, $
      /meanprofile, tracespec=1, traceorder=1, s2, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric
; ###########################################################################

; ###########################################################################
; second scan has higher S/N (28 versus 11)

    ispec, 'ngc_4100_drift_1_240.fits', /deredden, aperture=120, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_4100_drift_2_240.fits', /deredden, aperture=120, $
      /meanprofile, tracespec=1, traceorder=1, s2, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric
; ###########################################################################

; IUE galaxy; foreground star subtracted
    ispec, 'ngc_4102_drift_120.fits', /deredden, refrow=62, aperture=130.0, $
      /meanprofile, tracespec=1, traceorder=1, sbox=10, s1, $
      /starfit, star_refrow=93L, debug=debug, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; IUE galaxy
    ispec, 'ngc_4111_drift_180.fits', /deredden, aperture=110, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_4136_drift_180.fits', /deredden, aperture=200, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_4138_drift_060.fits', /deredden, aperture=120, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric
    
    ispec, 'ngc_4144_drift_180.fits', /deredden, aperture=110.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_4150_drift_060.fits', /deredden, aperture=100, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; dusty!    
    ispec, 'ngc_4157_drift_240.fits', /deredden, aperture=160.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_4163_drift_060.fits', /deredden, aperture=90.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; ###########################################################################
; the first scan has virtually no signal, so exclude from further
; analysis
    ispec, 'ngc_4183_drift_1_240.fits', /deredden, refrow=65, aperture=100, $
      /meanprofile, tracespec=1, traceorder=1, sbox=10, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=0, /noskymask, telluric=telluric

; foreground star subtracted
    ispec, 'ngc_4183_drift_2_240.fits', /deredden, refrow=65, aperture=100, $
      /meanprofile, tracespec=1, traceorder=1, sbox=10, s2, $
      /starfit, star_refrow=50L, frac_rows=0.04, debug=debug, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='ngc_4183_drift_240'
; ###########################################################################

    ispec, 'ngc_4190_drift_090.fits', /deredden, aperture=100.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; ###########################################################################
; GTO starburst

; Turner
    ispec, 'ngc_4194_drift_060.fits', /deredden, aperture=60.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_4194_drift_030.fits', /deredden, aperture=125.0, $
      /meanprofile, tracespec=1, traceorder=1, s2, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric
; ###########################################################################
    
; IUE galaxy; i think the bright point source is a star cluster within
; the galaxy, and not a foreground star
    ispec, 'ngc_4214_drift_120.fits', /deredden, refrow=60, aperture=190.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
;     /starfit, star_refrow=65L, frac_rows=0.04, trace_searchbox=4L, debug=debug, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; foreground star subtracted
    ispec, 'ngc_4217_drift_180.fits', /deredden, refrow=82, aperture=95.0, $
      /meanprofile, tracespec=1, traceorder=1, sbox=10, fluxrange=[-50,400]*1E-17, s1, $
      /starfit, star_refrow=100, frac_rows=0.05, debug=debug, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_4218_drift_040.fits', /deredden, aperture=70, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric
    
    ispec, 'ngc_4220_drift_120.fits', /deredden, aperture=150, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; foreground stellar contamination; unable to subtract the weak star     
    ispec, 'ngc_4244_drift_840.fits', /deredden, refrow=57, loaperture=60.0, upaperture=45.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; scan avoids star    
    ispec, 'ngc_4248_drift_090.fits', /deredden, loaperture=45.0, upaperture=35.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; foreground star subtracted
    ispec, 'ngc_4254_drift_240.fits', /deredden, refrow=80, loaperture=130.0, upaperture=110.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      /starfit, star_refrow=19L, frac_rows=0.05, debug=debug, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; NGC4258 --> IUE_EXTRACT    
    
    ispec, 'ngc_4288_drift_090.fits', /deredden, aperture=100, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; scan avoids star; foreground star excluded    
    ispec, 'ngc_4303_drift_240.fits', /deredden, refrow=65, loaperture=130.0, upaperture=65.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; M100; low-level foreground stellar contamination 
    ispec, 'ngc_4321_drift_360.fits', /deredden, aperture=180, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; foreground star subtracted    
    ispec, 'ngc_4384_drift_040.fits', /deredden, refrow=65, aperture=80.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      /starfit, star_refrow=45, frac_rows=0.05, debug=debug, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_4385_drift_060.fits', /deredden, aperture=100, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; ###########################################################################
; could use NGC4389 for aperture effects.  the second 75" scan has
; higher S/N (29 versus 21)
    ispec, 'ngc_4389_drift_060.fits', /deredden, aperture=140.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_4389_drift_1_075.fits', /deredden, aperture=140.0, $
      /meanprofile, tracespec=1, traceorder=1, s2, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_4389_drift_2_075.fits', /deredden, aperture=140.0, $
      /meanprofile, tracespec=1, traceorder=1, s3, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric
; ###########################################################################
    
    ispec, 'ngc_4414_drift_180.fits', /deredden, aperture=120.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_4455_drift_120.fits', /deredden, aperture=70.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; IUE galaxy
    ispec, 'ngc_4500_drift_080.fits', /deredden, refrow=64, aperture=80.0, $
      /meanprofile, tracespec=1, traceorder=1, fluxrange=[0.01,50.0]*1E-16, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; HOLM 419 B superposed on SE side
    ispec, 'ngc_4534_drift_060.fits', /deredden, aperture=120.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; IUE galaxy
    ispec, 'ngc_4569_drift_360.fits', /deredden, aperture=190.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_4605_drift_230.fits', /deredden, aperture=100.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_4618_drift_180.fits', /deredden, aperture=180, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric
    
    ispec, 'ngc_4625_drift_090.fits', /deredden, aperture=110, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric
    
    ispec, 'ngc_4651_drift_150.fits', /deredden, aperture=130, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_4656_drift_340.fits', /deredden, aperture=110, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_4666_drift_250.fits', /deredden, aperture=80, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; GTO starbust
    ispec, 'ngc_4670_drift_040.fits', /deredden, aperture=80.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; ###########################################################################
; Turner; two drift-scans were taken, one at 90 degrees with
; significant spatial overlap, and one at 180 degrees, which allows
; separation of the two individual galaxies; NGC4676 = ARP242

; POSANGLE=90; NGC4676A
    ispec, 'ngc_4676_drift_150.fits', /deredden, refrow=130, loaperture=8.0, upaperture=32.0, $
      /meanprofile, tracespec=1, traceorder=1, sbox=10, s1, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='ngc_4676a_drift_150', galaxy='NGC4676A'

; POSANGLE=90; NGC4676B    
    ispec, 'ngc_4676_drift_150.fits', /deredden, refrow=100, loaperture=30.0, upaperture=5.0, $
      /meanprofile, tracespec=1, traceorder=1, sbox=10, s2, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='ngc_4676b_drift_150', galaxy='NGC4676B'

; POSANGLE=90; NGC4676 (together)
    ispec, 'ngc_4676_drift_150.fits', /deredden, loaperture=43.0, upaperture=32.0, $
      /meanprofile, tracespec=1, traceorder=1, sbox=10, s3, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='ngc_4676_drift_150', galaxy='NGC4676'

; POSANGLE=180; NGC4676A
    ispec, 'ngc_4676_drift_155.fits', /deredden, refrow=150, loaperture=15.0, upaperture=30.0, $
      /meanprofile, tracespec=1, traceorder=1, sbox=10, s4, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='ngc_4676a_drift_155', galaxy='NGC4676A'

; POSANGLE=180; NGC4676B    
    ispec, 'ngc_4676_drift_155.fits', /deredden, refrow=110, loaperture=30.0, upaperture=20.0, $
      /meanprofile, tracespec=1, traceorder=1, sbox=10, s5, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='ngc_4676b_drift_155', galaxy='NGC4676B'

; POSANGLE=180; NGC4676 (together)
    ispec, 'ngc_4676_drift_155.fits', /deredden, loaperture=30.0, upaperture=65.0, $
      /meanprofile, tracespec=1, traceorder=1, s6, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='ngc_4676_drift_155', galaxy='NGC4676'
; ###########################################################################

    ispec, 'ngc_4713_drift_120.fits', /deredden, aperture=110, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; IUE galaxy    
    ispec, 'ngc_4736_drift_120.fits', /deredden, aperture=180, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; IUE galaxy; GTO starburst
    ispec, 'ngc_4861_drift_120.fits', /deredden, refrow=65, loaperture=50.0, upaperture=40.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; scan avoids star    
    ispec, 'ngc_4900_drift_060.fits', /deredden, aperture=110.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; ###########################################################################
; interacting pair; 180 degree position angle; extract both components
; and the N and S components separately; the N component is at smaller
; row number

; together
    ispec, 'ngc_4922_drift_040.fits', /deredden, loaperture=50.0, upaperture=30.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      galaxy='NGC4922'

; North; NGC4922 NED02
    ispec, 'ngc_4922_drift_040.fits', /deredden, refrow=53, loaperture=30.0, upaperture=5.0, $
      /meanprofile, tracespec=1, traceorder=1, sbox=5, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='ngc_4922n_drift_040', galaxy='NGC4922N'
    
; South; NGC4922 NED01
    ispec, 'ngc_4922_drift_040.fits', /deredden, refrow=68, loaperture=15.0, upaperture=30.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='ngc_4922s_drift_040', galaxy='NGC4922S'
; ###########################################################################

    ispec, 'ngc_5014_drift_040.fits', /deredden, aperture=90, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; ###########################################################################
; repeat observations

    ispec, 'ngc_5104_drift_1_060.fits', /deredden, aperture=40.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_5104_drift_2_060.fits', /deredden, aperture=40.0, $
      /meanprofile, tracespec=1, traceorder=1, s2, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric
; ###########################################################################

    ispec, 'ngc_5144_drift_060.fits', /deredden, loaperture=35.0, upaperture=30.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; M51; IUE galaxy; stitched together; foreground star subtracted
    ispec, 'ngc_5194_drift_360.fits', /deredden, aperture=375.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      /starfit, star_refrow=156, frac_rows=0.05, debug=debug, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; M83; crummy; exclude from atlas; foreground stellar contamination
    ispec, 'ngc_5236_drift_600.fits', /deredden, aperture=200.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=0, /noskymask, telluric=telluric

    ispec, 'ngc_5238_drift_050.fits', /deredden, aperture=85, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; IUE galaxy; GTO starburst
    ispec, 'ngc_5253_drift_080.fits', /deredden, aperture=130.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; IUE galaxy; interacting pair; unable to extract separate components
    ispec, 'ngc_5256_drift_030.fits', /deredden, aperture=50.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; ###########################################################################
; interacting pair; 90-degree slit position angle; NGC5257 is West of
; NGC5258

    ispec, 'ngc_5257-8_drift_090.fits', /deredden, refrow=89, loaperture=40.0, upaperture=30.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='ngc_5257_drift_090', galaxy='NGC5257'

    ispec, 'ngc_5257-8_drift_090.fits', /deredden, refrow=43, loaperture=35.0, upaperture=45.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='ngc_5258_drift_090', galaxy='NGC5258'

; together    
    ispec, 'ngc_5257-8_drift_090.fits', /deredden, refrow=89, loaperture=120.0, upaperture=30.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='arp_240_drift_090', galaxy='ARP240'
; ###########################################################################

; low-level foreground stellar contamination included 
    ispec, 'ngc_5264_drift_080.fits', /deredden, refrow=65, aperture=110.0, $
      /meanprofile, tracespec=1, traceorder=1, fluxrange=[-0.2,1.3]*1E-15, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; ###########################################################################
; NGC5278-9=ARP239; interacting pair; extract separately and together;
; significant spatial overlap 
    
; NGC5278
    ispec, 'ngc_5278-9_drift_042.fits', /deredden, refrow=155, loaperture=25.0, upaperture=30.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='ngc_5278_drift_042', galaxy='NGC5278'

; NGC5279
    ispec, 'ngc_5278-9_drift_042.fits', /deredden, refrow=110, loaperture=20.0, upaperture=12.0, $
      /meanprofile, tracespec=1, traceorder=1, sbox=10L, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='ngc_5279_drift_042', galaxy='NGC5279'

; ARP239 (together)
    ispec, 'ngc_5278-9_drift_042.fits', /deredden, loaperture=60.0, upaperture=30.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='arp_239_drift_042', galaxy='ARP239'
; ###########################################################################

    ispec, 'ngc_5430_drift_060.fits', /deredden, aperture=70.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric
    
; ###########################################################################
; interacting pair; 90 degree position angle; extract the sum
; (NGC5591), the East component (NGC5591E) and the West component
; (NGC5591W) separately; comparable S/N so choose the first scan 

; together
    ispec, 'ngc_5591_drift_1_030.fits', /deredden, refrow=58, loaperture=50.0, upaperture=30.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      galaxy='NGC5591'

; together
    ispec, 'ngc_5591_drift_2_030.fits', /deredden, refrow=58, loaperture=50.0, upaperture=30.0, $
      /meanprofile, tracespec=1, traceorder=1, s2, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      galaxy='NGC5591'

; East; NGC5591NED02
    ispec, 'ngc_5591_drift_1_030.fits', /deredden, refrow=44, loaperture=27.0, upaperture=7.0, $
      /meanprofile, tracespec=1, traceorder=1, sbox=10, s2, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='ngc_5591e_drift_030', galaxy='NGC5591E'
    
; West; NGC5591NED01
    ispec, 'ngc_5591_drift_1_030.fits', /deredden, refrow=68, loaperture=15.0, upaperture=30.0, $
      /meanprofile, tracespec=1, traceorder=1, sbox=10, s2, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='ngc_5591w_drift_030', galaxy='NGC5591W'
; ###########################################################################

    ispec, 'ngc_5607_drift_040.fits', /deredden, aperture=55.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_5653_drift_030.fits', /deredden, aperture=80.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric
    
; unable to subtract the (weak) foreground due to a poor trace at blue
; wavelengths 
    ispec, 'ngc_5676_drift_170.fits', /deredden, aperture=160.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
;     /starfit, star_refrow=82, frac_rows=0.05, trace_searchbox=2L, debug=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; ###########################################################################
; interacting pair, mildly blended; 90-degree slit position angle;
; NGC5953 (GTO starburst) is SW of NGC5954

; NGC5953
    ispec, 'ngc_5953-4_drift_075.fits', /deredden, refrow=85, loaperture=18.0, upaperture=42.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='ngc_5953_drift_075', galaxy='NGC5953'

; NGC5954    
    ispec, 'ngc_5953-4_drift_075.fits', /deredden, refrow=50, loaperture=22.0, upaperture=23.0, $
      /meanprofile, tracespec=1, traceorder=1, s2, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='ngc_5954_drift_075', galaxy='NGC5954'

; ARP091
    ispec, 'ngc_5953-4_drift_075.fits', /deredden, refrow=85, loaperture=63.0, upaperture=42.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='arp_091_drift_075', galaxy='ARP091'
; ###########################################################################

; foreground stars excluded    
    ispec, 'ngc_5992_drift_040.fits', /deredden, aperture=40.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; ###########################################################################
; interacting pair = ARP 072; NGC5996 is East (lower row number)
; relative to NGC5994; IUE galaxy; do not extract the pair because
; NGC5994 is of low S/N and the galaxies are well-separated spatially 

; NGC5996    
    ispec, 'ngc_5996_drift_120.fits', /deredden, loaperture=30.0, upaperture=50.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      galaxy='NGC5996'

; NGC5994; exclude from atlas -- low S/N
    ispec, 'ngc_5996_drift_120.fits', /deredden, refrow=115, aperture=30.0, $
      /meanprofile, tracespec=1, traceorder=1, fluxrange=[-5,100.0]*1E-17, s2, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=0, /noskymask, telluric=telluric, $
      outname='ngc_5994_drift_120', galaxy='NGC5994'
; ###########################################################################

; ###########################################################################
; IUE galaxy; spatially overlapping interacting pair; extract together 
    ispec, 'ngc_6052_drift_060.fits', /deredden, refrow=66, aperture=50.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric
; ###########################################################################

; ###########################################################################
; use NGC6090 to study aperture effects.  IUE galaxy; interacting
; pair; unable to extract separate components, so extract together 

    ispec, 'ngc_6090_drift_010.fits', /deredden, aperture=45.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_6090_drift_020.fits', /deredden, aperture=45.0, $
      /meanprofile, tracespec=1, traceorder=1, s2, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric
; ###########################################################################

    ispec, 'ngc_6207_drift_120.fits', /deredden, aperture=90, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; exclude foreground star    
    ispec, 'ngc_6240_drift_080.fits', /deredden, refrow=63, aperture=50.0, $ ; loaperture=25.0, upaperture=40.0, $
      /meanprofile, tracespec=1, traceorder=1, sbox=7, fluxrange=[-0.5,4]*1E-15, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; ###########################################################################
; ARP293 interacting pair; well-separated spatially, so do not extract
; the composite spectrum
    ispec, 'ngc_6285_drift_020.fits', /deredden, aperture=60, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric
 
   ispec, 'ngc_6286_drift_060.fits', /deredden, aperture=80, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric
; ###########################################################################

; ###########################################################################
; interacting triple (NGC6670A is East of NGC6670B)

; NGC6670A
    ispec, 'ngc_6670_drift_030.fits', /deredden, refrow=53, loaperture=25.0, upaperture=15.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='ngc_6670a_drift_030', galaxy='NGC6670A'

; NGC6670B
    ispec, 'ngc_6670_drift_030.fits', /deredden, refrow=70, loaperture=18.0, upaperture=22.0, $
      /meanprofile, tracespec=1, traceorder=1, sbox=10, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='ngc_6670b_drift_030', galaxy='NGC6670B'

; NGC6670 (together)
    ispec, 'ngc_6670_drift_030.fits', /deredden, refrow=53, loaperture=25.0, upaperture=55.0, $
      /meanprofile, tracespec=1, traceorder=1, sbox=10, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='ngc_6670_drift_030', galaxy='NGC6670'
; ###########################################################################

; foreground star excluded    
    ispec, 'ngc_6701_drift_060.fits', /deredden, refrow=62, loaperture=30.0, upaperture=45.0, $
      /meanprofile, tracespec=1, traceorder=1, sbox=10, fluxrange=[-0.1,2.1]*1E-15, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric
    
; ###########################################################################
; repeat observations; asymmetric spatial profile; the second scan has
; higher S/N; do not trace

    ispec, 'ngc_6926_drift_1_060.fits', /deredden, refrow=61, aperture=115.0, $
      /meanprofile, tracespec=0, traceorder=1, s1, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_6926_drift_2_060.fits', /deredden, refrow=61, aperture=115.0, $
      /meanprofile, tracespec=0, traceorder=1, s2, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric
; ###########################################################################

; IUE galaxy
    ispec, 'ngc_7130_drift_040.fits', /deredden, aperture=80, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; foreground star subtracted    
    ispec, 'ngc_7137_drift_060.fits', /deredden, refrow=56, aperture=80.0, $ ; loaperture=42.0, upaperture=23.0, $
      /meanprofile, tracespec=1, traceorder=1, sbox=10, s1, $
      /starfit, star_refrow=75, frac_rows=0.04, trace_searchbox=4L, debug=debug, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_7244_drift_040.fits', /deredden, aperture=35.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; ###########################################################################
; NGC7284-5=ARP093; interacting pair; significant spatial overlap

; NGC7284
    ispec, 'ngc_7284-5_drift_070.fits', /deredden, refrow=130, loaperture=7.0, upaperture=43.0, $
      /meanprofile, tracespec=1, traceorder=1, sbox=10, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='ngc_7284_drift_070', galaxy='NGC7284'

; NGC7285    
    ispec, 'ngc_7284-5_drift_070.fits', /deredden, refrow=110, loaperture=44.0, upaperture=6.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='ngc_7285_drift_070', galaxy='NGC7285'

; ARP093 (together)      
    ispec, 'ngc_7284-5_drift_070.fits', /deredden, loaperture=45.0, upaperture=55.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='arp_093_drift_070', galaxy='ARP093'
; ###########################################################################

; foreground star subtracted
    ispec, 'ngc_7316_drift_040.fits', /deredden, aperture=65.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      /starfit, star_refrow=51, frac_rows=0.04, trace_searchbox=4L, debug=debug, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_7448_drift_120.fits', /deredden, aperture=75.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric
    
    ispec, 'ngc_7465_drift_060.fits', /deredden, aperture=50.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_7468_drift_030.fits', /deredden, aperture=45.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; Seyfert 2! foreground star excluded
    ispec, 'ngc_7469_drift_060.fits', /deredden, refrow=64, loaperture=30.0, upaperture=35.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_7518_drift_040.fits', /deredden, aperture=75, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_7580_drift_040.fits', /deredden, refrow=60, aperture=50.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; ARP223    
    ispec, 'ngc_7585_drift_060.fits', /deredden, aperture=130.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; scan avoids star    
    ispec, 'ngc_7591_drift_075.fits', /deredden, aperture=80.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; ###########################################################################
; interacting pair (triple according to LEDA); significant spatial
; overlap, but by constraining the trace we can pull out three spectra 

; NGC7592B = NGC7592 EAST
    ispec, 'ngc_7592_drift_060.fits', /deredden, refrow=55, loaperture=22.0, upaperture=5.0, $
      /meanprofile, tracespec=1, traceorder=1, sbox=5, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='ngc_7592b_drift_060', galaxy='NGC7592B'

; NGC7592A = NGC7592 WEST
    ispec, 'ngc_7592_drift_060.fits', /deredden, refrow=61, loaperture=6.0, upaperture=25.0, $
      /meanprofile, tracespec=1, traceorder=1, trace_mincol=550, sbox=4, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='ngc_7592a_drift_060', galaxy='NGC7592A'

; NGC7592 (together)
    ispec, 'ngc_7592_drift_060.fits', /deredden, refrow=55, loaperture=22.0, upaperture=36.0, $
      /meanprofile, tracespec=1, traceorder=1, sbox=5, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      galaxy='NGC7592'

; ###########################################################################

    ispec, 'ngc_7620_drift_040.fits', /deredden, aperture=60, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric
    
    ispec, 'ngc_7624_drift_050.fits', /deredden, refrow=65, aperture=65.0, $
      /meanprofile, tracespec=1, traceorder=1, fluxrange=[-0.2,2.0]*1E-15, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_7625_drift_060.fits', /deredden, aperture=80, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; foreground stars subtracted; possibly some low-level foreground
; contamination 
    ispec, 'ngc_7640_drift_270.fits', /deredden, refrow=60, aperture=110.0, $
      /meanprofile, tracespec=1, traceorder=1, sbox=10, s1, $
      /starfit, star_refrow=[36,77], frac_rows=0.04, debug=debug, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_7673_drift_040.fits', /deredden, aperture=70.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; ###########################################################################
; NGC7674/NGC7674A = ARP 182; interacting pair; NGC7674
; (=MCG+01-59-080), a Sy2, is interacting with NGC7674A
; (=MCG+01-59-081), a companion dwarf galaxy

; NGC7674; foreground stellar contamination that is too weak to
; subtract 
    ispec, 'ngc_7674_drift_060.fits', /deredden, loaperture=24.0, upaperture=36.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
;     /starfit, star_refrow=50, frac_rows=0.06, debug=1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      galaxy='NGC7674'

; NGC7674A
    ispec, 'ngc_7674_drift_060.fits', /deredden, refrow=45, loaperture=18.0, upaperture=7.0, $
      /meanprofile, tracespec=1, traceorder=1, sbox=5, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='ngc_7674a_drift_060', galaxy='NGC7674A'

; ARP 182 (together); foreground stellar contamination   
    ispec, 'ngc_7674_drift_060.fits', /deredden, loaperture=49.0, upaperture=36.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='arp_182_drift_060', galaxy='ARP182'
; ###########################################################################
    
    ispec, 'ngc_7677_drift_060.fits', /deredden, aperture=50.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_7678_drift_090.fits', /deredden, aperture=105.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric
    
    ispec, 'ngc_7679_drift_030.fits', /deredden, aperture=80.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; low-level foreground stellar contamination    
    ispec, 'ngc_7713_drift_210.fits', /deredden, refrow=65, aperture=120.0, $
      /meanprofile, tracespec=1, traceorder=1, fluxrange=[-100,400]*1E-17, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_7714_drift_060.fits', /deredden, aperture=60, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; ARP222    
    ispec, 'ngc_7727_drift_100.fits', /deredden, aperture=150.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_7742_drift_060.fits', /deredden, refrow=50, aperture=100.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; foreground stellar contamination; unable to get a good subtraction
; due to low counts in the blue (the trace is okay, but the Gaussian
; fits are poor)
    ispec, 'ngc_7771_drift_050.fits', /deredden, aperture=130, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
;     /starfit, star_refrow=45, frac_rows=0.06, debug=debug, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_7782_drift_090.fits', /deredden, aperture=75.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_7798_drift_060.fits', /deredden, aperture=85, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ngc_7800_drift_060.fits', /deredden, loaperture=40.0, upaperture=60.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ugc_00591_drift_030.fits', /deredden, aperture=40, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ugc_00685_drift_040.fits', /deredden, aperture=80.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; dusty!
    ispec, 'ugc_00903_drift_060.fits', /deredden, aperture=90, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; foreground stellar contamination; unable to subtract 
    ispec, 'ugc_01281_drift_240.fits', /deredden, aperture=50.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ugc_01385_drift_030.fits', /deredden, aperture=45.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ugc_01561_drift_040.fits', /deredden, loaperture=30.0, upaperture=40.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ugc_02238_drift_030.fits', /deredden, aperture=60.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; ###########################################################################
; interacting pair; 180-degree slit position angle (see October 1998
; logs); North is "down" on the CCD; blended

; UGC02369 NED01; SOUTH; MCG+02-08-029
    ispec, 'ugc_02369_drift_040.fits', /deredden, refrow=63, loaperture=10.0, upaperture=20.0, $
      /meanprofile, tracespec=1, traceorder=1, sbox=10, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='ugc_02369s_drift_040', galaxy='UGC02369S'

; UGC02369 NED02; NORTH; MCG+02-08-030
    ispec, 'ugc_02369_drift_040.fits', /deredden, refrow=46, loaperture=18.0, upaperture=12.0, $
      /meanprofile, tracespec=1, traceorder=1, sbox=10, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='ugc_02369n_drift_040', galaxy='UGC02369N'

; UGC02369 (together)
    ispec, 'ugc_02369_drift_040.fits', /deredden, refrow=46, loaperture=18.0, upaperture=42.0, $
      /meanprofile, tracespec=1, traceorder=1, sbox=10, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='ugc_02369_drift_040', galaxy='UGC02369'
; ###########################################################################

; dusty!
    ispec, 'ugc_02982_drift_040.fits', /deredden, refrow=55, aperture=65, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric
    
; IUE galaxy
    ispec, 'ugc_03838_drift_040.fits', /deredden, aperture=40, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; SINGS; DDO053
    ispec, 'ugc_04459_drift_055.fits', /deredden, loaperture=50.0, upaperture=30.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; GTO starburst
    ispec, 'ugc_04483_drift_040.fits', /deredden, refrow=65, aperture=50.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric
    
; ARP195; galaxy triple; extract together; H-alpha is redshifted
; beyond the spectral range, so exclude from the atlas
    ispec, 'ugc_04653_drift_084.fits', /deredden, loaperture=50.0, upaperture=34.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=0, /noskymask, telluric=telluric

    ispec, 'ugc_04787_drift_120.fits', /deredden, aperture=40.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; interacting multiple; "Grasshopper" Galaxy; ARP 55; 90 degree
; position angle; extract together
    ispec, 'ugc_04881_drift_040.fits', /deredden, aperture=60, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; ###########################################################################
; interacting pair; ARP 300; 90 degree position angle; UGC05029 is
; East of UGC05028; the drift scan does not include all of UGC05029

; UGC05028    
    ispec, 'ugc_05028_drift_030.fits', /deredden, aperture=40, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='ugc_05028_drift_030', galaxy='UGC05028'
; ###########################################################################
    
; foreground star subtracted
    ispec, 'ugc_05101_drift_030.fits', /deredden, refrow=65, aperture=50.0, $
      /meanprofile, tracespec=1, traceorder=1, sbox=5, s1, $
      /starfit, star_refrow=55, frac_rows=0.06, debug=debug, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ugc_05151_drift_030.fits', /deredden, aperture=50.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; also called M81DWB; SINGS
    ispec, 'ugc_05423_drift_055.fits', /deredden, aperture=55.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; IUE galaxy (MRK33); SINGS
    ispec, 'ugc_05720_drift_040.fits', /deredden, aperture=50.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; ###########################################################################
; interacting pair; extract together    
    ispec, 'ugc_05998_drift_030.fits', /deredden, aperture=75, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; UGC05998 NED01; MRK0156    
    ispec, 'ugc_05998_drift_030.fits', /deredden, refrow=70, loaperture=15.0, upaperture=35.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      galaxy='UGC05998W'

; UGC05998 NED02
    ispec, 'ugc_05998_drift_030.fits', /deredden, refrow=55, loaperture=20.0, upaperture=5.0, $
      /meanprofile, tracespec=1, traceorder=1, sbox=5, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      galaxy='UGC05998E'

; UGC05998 (together)
    ispec, 'ugc_05998_drift_030.fits', /deredden, refrow=70, loaperture=40.0, upaperture=35.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      galaxy='UGC05998'
; ###########################################################################

    ispec, 'ugc_06029_drift_030.fits', /deredden, aperture=50.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ugc_06103_drift_030.fits', /deredden, aperture=50, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ugc_06399_drift_120.fits', /deredden, aperture=90, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; low-level foreground stellar contamination    
    ispec, 'ugc_06446_drift_120.fits', /deredden, aperture=65.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; GTO starburst; MRK0170; foreground star excluded
    ispec, 'ugc_06448_drift_040.fits', /deredden, aperture=60.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; IUE galaxy; GTO starburst; VIIZW403 
    ispec, 'ugc_06456_drift_060.fits', /deredden, refrow=63, loaperture=20.0, upaperture=30.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ugc_06520_drift_040.fits', /deredden, aperture=45, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; interacting pair; extract together; MRK178
    ispec, 'ugc_06541_drift_030.fits', /deredden, loaperture=35.0, upaperture=30.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; foreground star subtracted; foreground stellar contamination
    ispec, 'ugc_06628_drift_120.fits', /deredden, refrow=63, aperture=170.0, $
      /meanprofile, tracespec=1, traceorder=1, fluxrange=[-5,60]*1E-17, s1, $
      /starfit, star_refrow=30, debug=debug, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; GTO starburst; also called UM448
    ispec, 'ugc_06665_drift_060.fits', /deredden, aperture=30.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; foreground star excluded
    ispec, 'ugc_06667_drift_180.fits', /deredden, aperture=35, $
      /meanprofile, tracespec=1, traceorder=1, fluxrange=[-10,160]*1E-17, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; foreground star excluded    
    ispec, 'ugc_06773_drift_040.fits', /deredden, refrow=63, loaperture=52.0, upaperture=43.0, $
      /meanprofile, tracespec=1, traceorder=1, fluxrange=[-4,40]*1E-17, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ugc_06816_drift_090.fits', /deredden, aperture=90.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ugc_06818_drift_030.fits', /deredden, loaperture=50.0, upaperture=70.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; GTO starburst; also called UM462
    ispec, 'ugc_06850_drift_040.fits', /deredden, aperture=36.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ugc_06894_drift_015.fits', /deredden, refrow=60, aperture=90, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; very bright stellar contamination (late-type star Ca H+K); not a
; great subtraction
    ispec, 'ugc_06917_drift_120.fits', /deredden, refrow=55, aperture=150, $
      /meanprofile, tracespec=1, traceorder=1, fluxrange=[-0.1,2]*1E-15, s1, $
      /starfit, star_refrow=60, debug=debug, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ugc_06922_drift_040.fits', /deredden, refrow=60, aperture=80.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ugc_06923_drift_040.fits', /deredden, aperture=120.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ugc_06930_drift_120.fits', /deredden, aperture=160, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; cosmic ray!
    ispec, 'ugc_06969_drift_040.fits', /deredden, aperture=100, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric
    
    ispec, 'ugc_06983_drift_120.fits', /deredden, refrow=65, aperture=130, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; low-level foreground stellar contamination
    ispec, 'ugc_07089_drift_150.fits', /deredden, refrow=50, loaperture=70.0, upaperture=60.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ugc_07094_drift_040.fits', /deredden, aperture=60.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ugc_07354_drift_030.fits', /deredden, aperture=45.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ugc_07690_drift_060.fits', /deredden, aperture=90.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ugc_07950_drift_040.fits', /deredden, loaperture=25.0, upaperture=35.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; MRK 231 - whoa!
    ispec, 'ugc_08058_drift_040.fits', /deredden, aperture=40.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; DDO165
    ispec, 'ugc_08201_drift_120.fits', /deredden, loaperture=60, upaperture=80.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; giant HII region - not a foreground star!  (cf the 2D spectrum) 
    ispec, 'ugc_08323_drift_045.fits', /deredden, refrow=65, loaperture=35.0, upaperture=30.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; ###########################################################################
; interacting pair; 90-degree slit position angle; ARP238

; UGC08335 NED02; MCG+10-19-057; SE
    ispec, 'ugc_08335_drift_050.fits', /deredden, refrow=52, loaperture=20.0, upaperture=15.0, $
      /meanprofile, tracespec=1, traceorder=1, sbox=10, fluxrange=[-0.1,0.8]*1E-15, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='ugc_08335se_drift_050', galaxy='UGC08335SE'

; UGC08335 NED01; MCG+10-19-056; NW
    ispec, 'ugc_08335_drift_050.fits', /deredden, refrow=70, loaperture=15.0, upaperture=30.0, $
      /meanprofile, tracespec=1, traceorder=1, sbox=10, fluxrange=[-0.1,0.8]*1E-15, s2, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='ugc_08335nw_drift_050', galaxy='UGC08335NW'

; UGC08335 = ARP 238 (together)
    ispec, 'ugc_08335_drift_050.fits', /deredden, refrow=52, loaperture=20.0, upaperture=60.0, $
      /meanprofile, tracespec=1, traceorder=1, sbox=10, fluxrange=[-0.1,0.8]*1E-15, s3, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='ugc_08335_drift_050', galaxy='UGC08335'
; ###########################################################################
    
    ispec, 'ugc_08508_drift_060.fits', /deredden, aperture=90.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; MRK0273
    ispec, 'ugc_08696_drift_015.fits', /deredden, loaperture=60.0, upaperture=25.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; foreground star excluded
    ispec, 'ugc_09081_drift_040.fits', /deredden, loaperture=25.0, upaperture=30.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; GTO starburst; also called DDO187
    ispec, 'ugc_09128_drift_060.fits', /deredden, refrow=120, aperture=100.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ugc_09240_drift_090.fits', /deredden, aperture=110, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; ###########################################################################
; Turner; the SE component is UGC09425_NED02; the NW component is
; UGC09425_NED01; the composite system is called UGC09425=ARP241;
; spatial overlap

    ispec, 'ugc_09425_drift_025.fits', /deredden, refrow=140, loaperture=8.0, upaperture=27.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='ugc_09425nw_drift_025', galaxy='UGC09425NW'
    
    ispec, 'ugc_09425_drift_025.fits', /deredden, refrow=115, loaperture=27.0, upaperture=8.0, $
      /meanprofile, tracespec=1, traceorder=1, sbox=10, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='ugc_09425se_drift_025', galaxy='UGC09425SE'
    
    ispec, 'ugc_09425_drift_025.fits', /deredden, loaperture=45.0, upaperture=25.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='ugc_09425_drift_025', galaxy='UGC09425'
; ###########################################################################
    
; IUE galaxy
    ispec, 'ugc_09560_drift_030.fits', /deredden, aperture=40, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; ###########################################################################
; UGC09618 = ARP 302; interacting pair; repeat observations; stellar
; contamination; 0-degree position angle; the second scan has slightly
; higher S/N; do not include the foreground star in the extraction
; apertures; MCG+04-35-019 is North (higher row number) of MCG+04-35-018

; SOUTH; UGC09618NED01; MCG+04-35-018; foreground star subtracted
    ispec, 'ugc_09618_drift_1_030.fits', /deredden, refrow=46, loaperture=32.0, upaperture=15.0, $
      /meanprofile, tracespec=1, traceorder=1, sbox=7, s1, $
      /starfit, star_refrow=60, frac_rows=0.06, ngauss_terms=6L, debug=debug, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='ugc_09618s_1_drift_030', galaxy='UGC09618S'

; NORTH; UGC09618NED02; MCG+04-35-019; foreground star subtracted
    ispec, 'ugc_09618_drift_1_030.fits', /deredden, refrow=73, loaperture=10.0, upaperture=40.0, $
      /meanprofile, tracespec=1, traceorder=1, sbox=5, s2, $
      /starfit, star_refrow=60, frac_rows=0.06, ngauss_terms=6L, debug=debug, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='ugc_09618n_1_drift_030', galaxy='UGC09618N'

; UGC09618 = ARP 302 (together); foreground star subtracted
    ispec, 'ugc_09618_drift_1_030.fits', /deredden, refrow=73, loaperture=70.0, upaperture=40.0, $
      /meanprofile, tracespec=1, traceorder=1, sbox=5, s3, $
      /starfit, star_refrow=60, frac_rows=0.06, ngauss_terms=6L, debug=debug, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='ugc_09618_1_drift_030', galaxy='UGC09618'

; ---    
    
; SOUTH; UGC09618NED01; MCG+04-35-018; foreground star subtracted
    ispec, 'ugc_09618_drift_2_030.fits', /deredden, refrow=42, loaperture=30.0, upaperture=17.0, $
      /meanprofile, tracespec=1, traceorder=1, sbox=7, s4, $
      /starfit, star_refrow=60, frac_rows=0.06, ngauss_terms=6L, debug=debug, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='ugc_09618s_2_drift_030', galaxy='UGC09618S'

; NORTH; UGC09618NED02; MCG+04-35-019; foreground star subtracted
    ispec, 'ugc_09618_drift_2_030.fits', /deredden, refrow=73, loaperture=10.0, upaperture=40.0, $
      /meanprofile, tracespec=1, traceorder=1, sbox=5, s5, $
      /starfit, star_refrow=60, frac_rows=0.06, ngauss_terms=6L, debug=debug, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='ugc_09618n_2_drift_030', galaxy='UGC09618N'

; UGC09618 = ARP 302 (together); foreground star subtracted
    ispec, 'ugc_09618_drift_2_030.fits', /deredden, refrow=73, loaperture=70.0, upaperture=40.0, $
      /meanprofile, tracespec=1, traceorder=1, sbox=5, s6, $
      /starfit, star_refrow=60, frac_rows=0.06, ngauss_terms=6L, debug=debug, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='ugc_09618_2_drift_030', galaxy='UGC09618'

; ###########################################################################

; ###########################################################################
; Turner; UGC11175=ARP081=NGC6621-2; NGC6621 is NW of NGC6622;
; significant spatial overlap

; NGC6621; *NB* include the tidal debris feature in this spectrum  
    ispec, 'ugc_11175_drift_060.fits', /deredden, refrow=120, loaperture=27.0, upaperture=58.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='ngc_6621_drift_060', galaxy='NGC6621'

; NGC6622
    ispec, 'ugc_11175_drift_060.fits', /deredden, refrow=70, loaperture=27.0, upaperture=8.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='ngc_6622_drift_060', galaxy='NGC6622'

; UGC11175 (together)
    ispec, 'ugc_11175_drift_060.fits', /deredden, refrow=120, loaperture=62.0, upaperture=58.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='ugc_11175_drift_060', galaxy='UGC11175'
; ###########################################################################

; ###########################################################################
; interacting pair; 90 degree position angle

; UGC11680NED01 (IIZW101, MRK0897); WEST
    ispec, 'ugc_11680_drift_030.fits', /deredden, aperture=60.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric, $
      outname='ugc_11680w_drift_030', galaxy='UGC11680W'

; UGC11680NED02 (IIZW102, MRK0897); EAST; exclude from the
; atlas (jm03jul13uofa) 
    ispec, 'ugc_11680_drift_030.fits', /deredden, refrow=20, aperture=45, $
      /meanprofile, tracespec=1, traceorder=1, fluxrange=[-5,20]*1E-17, s2, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=0, /noskymask, telluric=telluric, $
      outname='ugc_11680e_drift_030', galaxy='UGC11680E'
; ###########################################################################
    
    ispec, 'ugc_12150_drift_040.fits', /deredden, refrow=62, aperture=55.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric
    
; ###########################################################################
; repeat observations; use UGC12490 to study aperture effects 

    ispec, 'ugc_12490_drift_030.fits', /deredden, refrow=50, aperture=50.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ugc_12490_drift_040.fits', /deredden, refrow=50, aperture=50.0, $
      /meanprofile, tracespec=1, traceorder=1, s2, $
      datapath=datapath, outpath=repeatpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric
; ###########################################################################

; foreground star subtracted    
    ispec, 'ugc_12588_drift_045.fits', /deredden, refrow=65, aperture=100, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      /starfit, star_refrow=60, frac_rows=0.06, debug=debug, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ugc_12747_drift_040.fits', /deredden, loaperture=30.0, upaperture=25.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ugca_073_drift_060.fits', /deredden, aperture=40.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; foreground star subtracted; the star near row 62 is difficult to
; subtract so leave it in; stellar contamination
    ispec, 'ugca_090_drift_300.fits', /deredden, refrow=60, aperture=110, $
      /meanprofile, tracespec=1, traceorder=1, sbox=10, fluxrange=[-50,220]*1E-17, s1, $
      /starfit, star_refrow=80, frac_rows=0.06, debug=debug, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; foreground stellar contamination
    ispec, 'ugca_114_drift_120.fits', /deredden, refrow=62, aperture=100, $
      /meanprofile, tracespec=1, traceorder=1, fluxrange=[-15,140]*1E-17, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; GTO starburst; also alled II Zw 040
    ispec, 'ugca_116_drift_030.fits', /deredden, refrow=65, aperture=35, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; foreground star excluded
    ispec, 'ugca_130_drift_040.fits', /deredden, aperture=25.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; DDO44; S/N = 1!; remove from atlas (jm02nov1uofa)
    ispec, 'ugca_133_drift_120.fits', /deredden, refrow=60, aperture=60, $
      /meanprofile, tracespec=1, traceorder=1, fluxrange=[-5,50]*1E-17, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=0, /noskymask, telluric=telluric
    
; IUE galaxy; I Zw 018; interacting pair; extract together; GTO
; starburst
    ispec, 'ugca_166_drift_020.fits', /deredden, refrow=65, aperture=20.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ugca_208_drift_020.fits', /deredden, aperture=30.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; IUE galaxy; GTO starburst; also called MRK0153
    ispec, 'ugca_219_drift_030.fits', /deredden, aperture=35, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; IUE galaxy
    ispec, 'ugca_225_drift_020.fits', /deredden, aperture=25.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; IUE galaxy
    ispec, 'ugca_281_drift_030.fits', deredden=1, refrow=58, loaperture=35.0, upaperture=15.0, $
      /meanprofile, tracespec=1, traceorder=1, fluxrange=[-5,110]*1E-17, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric
    
; GTO starburst; S/N = 1; H-alpha marginally detected; exclude from
; atlas: low S/N (jm05jul27uofa)
    ispec, 'ugca_292_drift_060.fits', /deredden, refrow=135, aperture=50.0, $
      /meanprofile, tracespec=0, traceorder=1, fluxrange=[-5,30]*1E-17, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=0, /noskymask, telluric=telluric

    ispec, 'ugca_320_drift_090.fits', /deredden, aperture=140.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ugca_372_drift_015.fits', /deredden, aperture=20.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric
    
; IUE galaxy
    ispec, 'ugca_410_drift_015.fits', /deredden, aperture=25.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

    ispec, 'ugca_439_drift_020.fits', /deredden, aperture=25.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

; GTO starburst
    ispec, 'um_461_drift_020.fits', /deredden, aperture=25.0, $
      /meanprofile, tracespec=1, traceorder=1, s1, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /noskymask, telluric=telluric

;; combine repeaters - this has to be done after running ATLAS1D_NUCLEAR_EXTRACT!
;
;   splog, 'Combining repeaters.'
;   atlas1d_combine_repeaters, wfits=wfits
    
    splog, 'Total time to extract all spectra (minutes):', (systime(1)-t0)/60.0

; atlas1d_header_redshifts, /update    
; atlas1d_info, info, /write

return
end
