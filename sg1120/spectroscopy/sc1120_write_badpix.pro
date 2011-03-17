pro make_badpix
; jm04nov07uofa
; convert the FITS file bad pixel masks provided by Anthony into
; iSPEC2d-compatible text files

    rootpath = '/home/ioannis/research/projects/sc1120/'

    cwdpath = rootpath+'redux/'
    datapath = rootpath+'data/'

    q1 = transpose(readfits(cwdpath+'Q1.bpm.fits',h1,/silent))
;   q2 = transpose(readfits(cwdpath+'Q2.bpm.fits',h1,/silent))
;   q3 = transpose(readfits(cwdpath+'Q3.bpm.fits',h1,/silent))
;   q4 = transpose(readfits(cwdpath+'Q4.bpm.fits',h1,/silent))

    ncols = 4096L
    nrows = 2140L
    
    badpixfile = ['Q1_badpix.dat','Q2_badpix.dat','Q3_badpix.dat','Q4_badpix.dat']

    for j = 0L, 3L do begin
    
       openw, lun, badpixfile[j], /get_lun
       
       for i = 0L, ncols-1L do begin

          bad = where(q1[i,*] eq 1L,nbad)
                    
          
       endfor

       free_lun, lun

    endfor
       
return
end
    
