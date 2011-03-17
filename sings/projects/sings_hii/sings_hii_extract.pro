pro sings_hii_extract, noplot=noplot, wfits=wfits
; jm06jul31uofa - 

    m = mrdfits('/home/ioannis/kennicutt/sings/projects/log12oh/sings_hiiregions_v3.0.fits.gz',1)
    struct_print, struct_trimtags(m,select=['hii_galaxy','galaxy_ra','galaxy_dec','galaxy_rc3_r25',$
      'galaxy_rc3_pa','galaxy_rc3_incl','hii_region','hii_raoffset','hii_deoffset',$
      'hii_rc3_radius','hii_rc3_rr25','zstrong_r23','reference']), file='sings_hiiregions.dat'
;   mpage -1 -l -W105 sings_hiiregions.dat > dum

    datapath = sings_hii_path(/spec2d)
    outpath = sings_hii_path(/spec1d)

; ###########################################################################    
; NGC2403 - BRENT
; ###########################################################################    
    
; HK065
    ispec, 'ngc_2403_hk065.fits', /deredden, refrow=65, sbox=10, aperture=20.0, $
      /meanprofile, tracespec=1, traceorder=2, skyaperture=40.0, flanking=25.0, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask

; HK128
    ispec, 'ngc_2403_hk128.fits', /deredden, refrow=65, sbox=10, aperture=20.0, $
      /meanprofile, tracespec=1, traceorder=2, skyaperture=40.0, /flanking, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask

; HK270 - funny looking; clouds; maybe don't combine
;                 sequential images
    ispec, 'ngc_2403_hk270.fits', /deredden, refrow=65, sbox=10, aperture=15.0, $
      /meanprofile, tracespec=1, traceorder=2, skyaperture=30.0, /flanking, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask

; HK361    
    ispec, 'ngc_2403_hk361.fits', /deredden, refrow=65, sbox=10, aperture=20.0, $
      /meanprofile, tracespec=1, traceorder=2, skyaperture=40.0, /flanking, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask

; HK542
    ispec, 'ngc_2403_hk542.fits', /deredden, refrow=65, sbox=10, aperture=20.0, $
      /meanprofile, tracespec=1, traceorder=2, skyaperture=20.0, /flanking, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask

; VS9
    ispec, 'ngc_2403_vs9.fits', /deredden, refrow=65, sbox=10, aperture=15.0, $
      /meanprofile, tracespec=1, traceorder=2, skyaperture=30.0, /flanking, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask

; ###########################################################################    
; NGC3521 - BRENT
; ###########################################################################    
    
; P02
    ispec, 'ngc_3521_p02.fits', /deredden, refrow=59, sbox=10, aperture=10.0, $
      /meanprofile, tracespec=0, traceorder=1, skyaperture=20.0, /flanking, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask

; P05 - astrometry needed
    ispec, 'ngc_3521_p05.fits', /deredden, refrow=59, sbox=10, aperture=10.0, $
      /meanprofile, tracespec=0, traceorder=1, skyaperture=40.0, skylower=150.0, skyupper=0.0, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask

; P58 - check 2D reductions
    ispec, 'ngc_3521_p58.fits', /deredden, refrow=47, sbox=10, aperture=10.0, $
      /meanprofile, tracespec=0, traceorder=1, skyaperture=10.0, /flanking, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask

; ###########################################################################    
; NGC3627 - BRENT
; ###########################################################################    
    
; P04 - sky subtraction?  send to Rob?
    ispec, 'ngc_3627_p04.fits', /deredden, refrow=57, sbox=10, aperture=10.0, $
      /meanprofile, tracespec=0, traceorder=1, skyaperture=20.0, skylower=0.0, skyupper=0.0, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask

; P07 - compare 24-micron and H-alpha profile
    ispec, 'ngc_3627_p07.fits', /deredden, refrow=57, sbox=10, aperture=10.0, $
      /meanprofile, tracespec=0, traceorder=1, skyaperture=20.0, skylower=0.0, skyupper=0.0, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask

; P23 - compare 24-micron and H-alpha profile
    ispec, 'ngc_3627_p23.fits', /deredden, refrow=55, sbox=5, aperture=10.0, $
      /meanprofile, tracespec=0, traceorder=1, skyaperture=40.0, skylower=0.0, skyupper=0.0, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask

; ###########################################################################    
; NGC3938
; ###########################################################################    

; P27 - ?
    ispec, 'ngc_3938_p27.fits', /deredden, refrow=70, sbox=10, aperture=10.0, $
      /meanprofile, tracespec=0, traceorder=1, skyaperture=20.0, skylower=5.0, skyupper=5.0, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask

; P51
    ispec, 'ngc_3938_p51_p59.fits', /deredden, refrow=68, sbox=5, aperture=10.0, $
      /meanprofile, tracespec=1, traceorder=1, skyaperture=[50.0,20.0], skylower=5.0, skyupper=5.0, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, $
      outname='ngc_3938_p51.fits', galaxy='NGC3938/P51'      

; P59
    ispec, 'ngc_3938_p51_p59.fits', /deredden, refrow=78, sbox=5, aperture=10.0, $
      /meanprofile, tracespec=1, traceorder=1, skyaperture=[50.0,20.0], skylower=5.0, skyupper=5.0, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, $
      outname='ngc_3938_p59.fits', galaxy='NGC3938/P59'      

; ###########################################################################    
; NGC4254
; ###########################################################################    

; P02 - ?
    ispec, 'ngc_4254_p02_p08.fits', /deredden, refrow=70, sbox=10, aperture=10.0, $
      /meanprofile, tracespec=0, traceorder=1, skyaperture=20.0, skylower=5.0, skyupper=5.0, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, $
      outname='ngc_4254_p02.fits', galaxy='NGC4254/P02'

; P08 - ?
    ispec, 'ngc_4254_p02_p08.fits', /deredden, refrow=115, sbox=5, aperture=10.0, $
      /meanprofile, tracespec=0, traceorder=1, skyaperture=20.0, skylower=5.0, skyupper=5.0, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, $
      outname='ngc_4254_p08.fits', galaxy='NGC4254/P08'

; ###########################################################################    
; NGC4321
; ###########################################################################    

; P39

    ispec, 'ngc_4321_p39_1.fits', /deredden, refrow=57, sbox=10, aperture=10.0, $
      /meanprofile, tracespec=1, traceorder=1, skyaperture=20.0, /flanking, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask
    ispec, 'ngc_4321_p39_2.fits', /deredden, refrow=57, sbox=10, aperture=10.0, $
      /meanprofile, tracespec=1, traceorder=1, skyaperture=20.0, /flanking, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask

; P44

    ispec, 'ngc_4321_p44_1.fits', /deredden, refrow=57, sbox=10, aperture=10.0, $
      /meanprofile, tracespec=1, traceorder=1, skyaperture=20.0, /flanking, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask
    ispec, 'ngc_4321_p44_2.fits', /deredden, refrow=57, sbox=10, aperture=10.0, $
      /meanprofile, tracespec=1, traceorder=1, skyaperture=20.0, /flanking, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask

; P57    
    
    ispec, 'ngc_4321_p57.fits', /deredden, refrow=57, sbox=10, aperture=10.0, $
      /meanprofile, tracespec=1, traceorder=1, skyaperture=20.0, /flanking, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask    
    
; ###########################################################################    
; NGC4631 - BRENT
; ###########################################################################    

; P01
    
    ispec, 'ngc_4631_p01.fits', /deredden, refrow=68, sbox=10, aperture=10.0, $
      /meanprofile, tracespec=1, traceorder=1, skyaperture=40.0, skylower=0.0, skyupper=0.0, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask    
        
; P26
    
    ispec, 'ngc_4631_p26.fits', /deredden, refrow=68, sbox=10, aperture=10.0, $
      /meanprofile, tracespec=1, traceorder=1, skyaperture=40.0, skylower=0.0, skyupper=0.0, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask    
        
; P34
    
    ispec, 'ngc_4631_p34.fits', /deredden, refrow=68, sbox=10, aperture=10.0, $
      /meanprofile, tracespec=1, traceorder=1, skyaperture=[20.0,50.0], skylower=0.0, skyupper=0.0, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask    
        
; ###########################################################################    
; NGC5055 - BRENT
; ###########################################################################    

; SL01

    ispec, 'ngc_5055_sl01.fits', /deredden, refrow=70, sbox=10, aperture=10.0, $
      /meanprofile, tracespec=1, traceorder=1, skyaperture=20.0, /flanking, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask

; ###########################################################################    
; NGC5194 - BRENT
; ###########################################################################    

; CCM001 - ask Moire?

    ispec, 'ngc_5194_ccm001_ccm010.fits', /deredden, refrow=59, sbox=10, aperture=10.0, $
      /meanprofile, tracespec=1, traceorder=1, skyaperture=20.0, skylower=0.0, skyupper=0.0, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, $
      outname='ngc_5194_ccm001.fits', galaxy='NGC5194/CCM001'      

; CCM010 - ask Moire?

    ispec, 'ngc_5194_ccm001_ccm010.fits', /deredden, refrow=68, sbox=10, aperture=10.0, $
      /meanprofile, tracespec=1, traceorder=1, skyaperture=20.0, /flanking, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask, $
      outname='ngc_5194_ccm010.fits', galaxy='NGC5194/CCM010'      

; SL01 - compare 24-micron and H-alpha profile

    ispec, 'ngc_5194_sl01.fits', /deredden, refrow=57, sbox=10, aperture=10.0, $
      /meanprofile, tracespec=1, traceorder=1, skyaperture=20.0, skylower=0.0, skyupper=180.0, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask

; SL04 - ?

    ispec, 'ngc_5194_sl04.fits', /deredden, refrow=57, sbox=10, aperture=10.0, $
      /meanprofile, tracespec=0, traceorder=1, skyaperture=20.0, /flanking, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask

; ###########################################################################    
; NGC5713
; ###########################################################################    

; P04 - ?

    ispec, 'ngc_5713_p04_1.fits', /deredden, refrow=68, sbox=10, aperture=10.0, $
      /meanprofile, tracespec=1, traceorder=1, skyaperture=20.0, /flanking, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask
    ispec, 'ngc_5713_p04_2.fits', /deredden, refrow=68, sbox=10, aperture=10.0, $
      /meanprofile, tracespec=1, traceorder=1, skyaperture=20.0, /flanking, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask

; P05 - ?
    
    ispec, 'ngc_5713_p05.fits', /deredden, refrow=68, sbox=10, aperture=10.0, $
      /meanprofile, tracespec=1, traceorder=1, skyaperture=20.0, /flanking, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask

; ###########################################################################    
; NGC6946 - BRENT
; ###########################################################################    

; HK003
    
    ispec, 'ngc_6946_hk003.fits', /deredden, refrow=68, sbox=10, aperture=10.0, $
      /meanprofile, tracespec=1, traceorder=1, skyaperture=20.0, /flanking, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask

; HK004
    
    ispec, 'ngc_6946_hk004.fits', /deredden, refrow=68, sbox=10, aperture=10.0, $
      /meanprofile, tracespec=1, traceorder=1, skyaperture=20.0, /flanking, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask

; HK028
    
    ispec, 'ngc_6946_hk028.fits', /deredden, refrow=68, sbox=10, aperture=10.0, $
      /meanprofile, tracespec=1, traceorder=1, skyaperture=20.0, /flanking, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask

; HK040
    
    ispec, 'ngc_6946_hk040.fits', /deredden, refrow=68, sbox=10, aperture=10.0, $
      /meanprofile, tracespec=1, traceorder=1, skyaperture=20.0, /flanking, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask

; HK288
    
    ispec, 'ngc_6946_hk288.fits', /deredden, refrow=68, sbox=10, aperture=10.0, $
      /meanprofile, tracespec=1, traceorder=1, skyaperture=20.0, /flanking, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask

; SL01 - ?
    
    ispec, 'ngc_6946_sl01.fits', /deredden, refrow=59, sbox=10, aperture=10.0, $
      /meanprofile, tracespec=0, traceorder=1, skyaperture=20.0, /flanking, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask

; SL02 - ?
    
    ispec, 'ngc_6946_sl02.fits', /deredden, refrow=59, sbox=10, aperture=10.0, $
      /meanprofile, tracespec=0, traceorder=1, skyaperture=20.0, /flanking, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask

; SL03 - ?
    
    ispec, 'ngc_6946_sl03.fits', /deredden, refrow=59, sbox=10, aperture=10.0, $
      /meanprofile, tracespec=0, traceorder=1, skyaperture=[30.0,50.0], skylower=0.0, skyupper=0.0, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask

; SL04 - ?
    
    ispec, 'ngc_6946_sl04.fits', /deredden, refrow=59, sbox=10, aperture=10.0, $
      /meanprofile, tracespec=0, traceorder=1, skyaperture=20.0, /flanking, $
      datapath=datapath, outpath=outpath, noplot=noplot, wfits=wfits, /nobadmask, /noskymask

; ###########################################################################    

return
end
