pro delete_2003_chip2
; jm04jun20uofa
; remove chip2 files from the RAW directory

; also remove calibration data that does not apply to this project;
; use INS_SLIT_NAME tag.  the following objects were taken with the
; 0.3" slit

; FORS2.2003-03-27T10:51:49.329.fits
; FORS2.2003-03-27T10:53:39.199.fits
; FORS2.2003-03-27T10:54:58.540.fits
; FORS2.2003-03-27T10:56:16.174.fits
; FORS2.2003-03-27T10:57:36.049.fits
; FORS2.2003-03-27T10:58:55.915.fits
; FORS2.2003-03-28T10:40:32.695.fits
; FORS2.2003-03-28T10:42:27.920.fits
; FORS2.2003-03-28T10:43:46.455.fits
; FORS2.2003-03-28T10:45:04.997.fits
; FORS2.2003-03-28T10:46:22.927.fits
; FORS2.2003-03-28T10:47:39.897.fits
    
; verify the INS_GRIS1_NAME tag to verify the grating setting.  remove
; the following dome flats that were taken with the GRIS_150I grating: 

; FORS2.2003-03-27T11:00:22.489.fits
; FORS2.2003-03-27T11:01:19.072.fits
; FORS2.2003-03-27T11:02:16.640.fits
; FORS2.2003-03-27T11:03:13.778.fits
; FORS2.2003-03-27T12:42:05.634.fits
; FORS2.2003-03-28T10:49:06.887.fits
; FORS2.2003-03-28T10:50:01.037.fits
; FORS2.2003-03-28T10:50:55.556.fits
; FORS2.2003-03-28T10:51:49.981.fits
; FORS2.2003-03-28T12:34:53.837.fits
; FORS2.2003-03-29T13:13:05.186.fits
; FORS2.2003-03-30T12:23:49.682.fits
; FORS2.2003-03-31T13:12:44.467.fits
; FORS2.2003-04-01T13:11:10.923.fits
; FORS2.2003-04-02T13:01:23.764.fits
; FORS2.2003-04-03T13:02:47.163.fits
; FORS2.2003-04-04T12:37:01.876.fits
; FORS2.2003-04-05T12:48:31.284.fits
; FORS2.2003-04-06T13:16:10.216.fits
    
    rawpath = '/home/ioannis/ediscs/2003/raw/'
    
    flist = file_search(rawpath+'*.fits',count=fcount)

    for i = 0L, fcount-1L do begin

       h = headfits(flist[i])
       origfile = sxpar(h,'ORIGFILE')

       if strmatch(origfile,'*chip2*',/fold) then spawn, ['/bin/rm '+flist[i]], /sh ; print, i, origfile
    
    endfor

return
end
    
