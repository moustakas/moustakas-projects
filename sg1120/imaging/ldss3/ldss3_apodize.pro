function ldss3_apodize, imsize, file=file, chip1=chip1
; jm07jan09nyu - apodize LDSS3 images based on an elliptical aperture
;   that I defined "by hand" using DS9; this apodization image should
;   be applied to the inverse variance map *and* the data because it
;   also crops nicely

    if (n_elements(imsize) ne 2L) then begin
       splog, 'IMSIZE must be a two-element array.'
       return, 1.0
    endif

;   rmajor = 1250.0             ; 1285.0
;   rminor = 1240.0             ; 1253.0
;   e_pa = 0.0
;   e_ycen = imsize[1]/2.0
;
;   if keyword_set(chip1) then begin
;      e_xcen = 0.0
;   endif else begin
;      e_xcen = imsize[0]
;   endelse
;   
;   dist_ellipse, distimage, imsize, e_xcen, e_ycen, rmajor/rminor, 90.0+e_pa ; note angle offset wrt DS9

    if keyword_set(chip1) then begin
       rmajor = 1250.0
       e_xcen = 0.0
       e_ycen = imsize[1]/2.0 - 10.0 ; NOTE! Trim-dependent! 
    endif else begin
       rmajor = 1280.0
       e_xcen = imsize[0]
       e_ycen = imsize[1]/2.0 - 10.0 ; NOTE! Trim-dependent! 
    endelse
    
    dist_circle, distimage, imsize, e_xcen, e_ycen

; ---------------------------------------------------------------------------
; arctangent apodization    
; ---------------------------------------------------------------------------
    param = 50.0
;;  apimage = cos((!pi/2.0)*((distimage/rmajor) < 1.0))
    apimage = ((atan(param*(1.0-(distimage/rmajor)))*2.0/!pi) > 0.0) < 1.0
    apimage = temporary(apimage) * (distimage lt rmajor)
; ---------------------------------------------------------------------------

    apimage = (distimage lt rmajor)

; customized bad pixel masks; depends on the trim!!

; ---------------------------------------------------------------------------    
    if (n_elements(file) ne 0L) then if strmatch(file,'*2016_sg1120_2_g_c1*') then begin
       rmajor = 430.0 & rminor = 200.0 & e_pa = 0.0 & e_xcen = 1215.0 & e_ycen = 1345.0
       dist_ellipse, mask, imsize, e_xcen, e_ycen, rmajor/rminor, e_pa
       apimage = temporary(apimage) * (mask gt rmajor)
    endif
    
    if (n_elements(file) ne 0L) then if strmatch(file,'*2017_sg1120_2_r_c1*') then begin
       rmajor = 430.0 & rminor = 200.0 & e_pa = 0.0 & e_xcen = 1215.0 & e_ycen = 1345.0
       dist_ellipse, mask, imsize, e_xcen, e_ycen, rmajor/rminor, e_pa
       apimage = temporary(apimage) * (mask gt rmajor)
    endif

; ---------------------------------------------------------------------------    
    if (n_elements(file) ne 0L) then if strmatch(file,'*2018_sg1120_2_r_c1*') then begin
       rmajor = 492.0 & rminor = 270.0 & e_pa = 25.0 & e_xcen = 1095.0 & e_ycen = 1750.0
       dist_ellipse, mask, imsize, e_xcen, e_ycen, rmajor/rminor, e_pa
       apimage = temporary(apimage) * (mask gt rmajor)
    endif

    if (n_elements(file) ne 0L) then if strmatch(file,'*2019_sg1120_2_g_c1*') then begin
       rmajor = 492.0 & rminor = 270.0 & e_pa = 25.0 & e_xcen = 1095.0 & e_ycen = 1750.0
       dist_ellipse, mask, imsize, e_xcen, e_ycen, rmajor/rminor, e_pa
       apimage = temporary(apimage) * (mask gt rmajor)
    endif
    
; ---------------------------------------------------------------------------    
    if (n_elements(file) ne 0L) then if strmatch(file,'*2018_sg1120_2_r_c2*') then begin
       rmajor = 545.0 & rminor = 190.0 & e_pa = 330.0 & e_xcen = 930.0 & e_ycen = 375.0
       dist_ellipse, mask, imsize, e_xcen, e_ycen, rmajor/rminor, e_pa-90.0 ; note angle offset
       apimage = temporary(apimage) * (mask gt rmajor)
    endif

    if (n_elements(file) ne 0L) then if strmatch(file,'*2019_sg1120_2_g_c2*') then begin
       rmajor = 545.0 & rminor = 190.0 & e_pa = 330.0 & e_xcen = 930.0 & e_ycen = 375.0
       dist_ellipse, mask, imsize, e_xcen, e_ycen, rmajor/rminor, e_pa-90.0 ; note angle offset
       apimage = temporary(apimage) * (mask gt rmajor)
    endif

; ---------------------------------------------------------------------------    
    if (n_elements(file) ne 0L) then if strmatch(file,'*2020_sg1120_2_r_c1*') then begin
       rmajor = 820.0 & rminor = 525.0 & e_pa = 0.0 & e_xcen = 885.0 & e_ycen = 1235.0
       dist_ellipse, mask, imsize, e_xcen, e_ycen, rmajor/rminor, e_pa
       apimage = temporary(apimage) * (mask gt rmajor)
    endif

    if (n_elements(file) ne 0L) then if strmatch(file,'*2021_sg1120_2_g_c1*') then begin
       rmajor = 820.0 & rminor = 525.0 & e_pa = 0.0 & e_xcen = 885.0 & e_ycen = 1235.0
       dist_ellipse, mask, imsize, e_xcen, e_ycen, rmajor/rminor, e_pa
       apimage = temporary(apimage) * (mask gt rmajor)
    endif

; ---------------------------------------------------------------------------    
    if (n_elements(file) ne 0L) then if strmatch(file,'*2022_sg1120_3_g_c2*') then begin
       rmajor = 260.0 & rminor = 52.0 & e_pa = 35.0 & e_xcen = 465.0 & e_ycen = 640.0
       dist_ellipse, mask, imsize, e_xcen, e_ycen, rmajor/rminor, e_pa
       apimage = temporary(apimage) * (mask gt rmajor)
    endif

    if (n_elements(file) ne 0L) then if strmatch(file,'*2023_sg1120_3_r_c2*') then begin
       rmajor = 260.0 & rminor = 52.0 & e_pa = 35.0 & e_xcen = 465.0 & e_ycen = 640.0
       dist_ellipse, mask, imsize, e_xcen, e_ycen, rmajor/rminor, e_pa
       apimage = temporary(apimage) * (mask gt rmajor)
    endif

    if (n_elements(file) ne 0L) then if strmatch(file,'*2024_sg1120_3_r_c2*') then begin
       rmajor = 260.0 & rminor = 52.0 & e_pa = 35.0 & e_xcen = 465.0 & e_ycen = 640.0
       dist_ellipse, mask, imsize, e_xcen, e_ycen, rmajor/rminor, e_pa
       apimage = temporary(apimage) * (mask gt rmajor)
    endif

    if (n_elements(file) ne 0L) then if strmatch(file,'*2025_sg1120_3_g_c2*') then begin
       rmajor = 260.0 & rminor = 52.0 & e_pa = 35.0 & e_xcen = 465.0 & e_ycen = 640.0
       dist_ellipse, mask, imsize, e_xcen, e_ycen, rmajor/rminor, e_pa
       apimage = temporary(apimage) * (mask gt rmajor)
    endif

; ---------------------------------------------------------------------------    
    if (n_elements(file) ne 0L) then if strmatch(file,'*2032_sg1120_4_g_c2*') then begin
       rmajor = 550.0 & rminor = 320.0 & e_pa = 0.0 & e_xcen = 460.0 & e_ycen = 1670.0
       dist_ellipse, mask, imsize, e_xcen, e_ycen, rmajor/rminor, e_pa
       apimage = temporary(apimage) * (mask gt rmajor)
    endif

    if (n_elements(file) ne 0L) then if strmatch(file,'*2033_sg1120_4_r_c2*') then begin
       rmajor = 550.0 & rminor = 320.0 & e_pa = 0.0 & e_xcen = 460.0 & e_ycen = 1670.0
       dist_ellipse, mask, imsize, e_xcen, e_ycen, rmajor/rminor, e_pa
       apimage = temporary(apimage) * (mask gt rmajor)
    endif

; ---------------------------------------------------------------------------    
    if (n_elements(file) ne 0L) then if strmatch(file,'*2034_sdss_1_r_c1*') then begin
       rmajor = 500.0 & rminor = 240.0 & e_pa = 40.0 & e_xcen = 875.0 & e_ycen = 2075.0
       dist_ellipse, mask, imsize, e_xcen, e_ycen, rmajor/rminor, e_pa
       apimage = temporary(apimage) * (mask gt rmajor)
    endif

    if (n_elements(file) ne 0L) then if strmatch(file,'*2035_sdss_1_g_c1*') then begin
       rmajor = 500.0 & rminor = 240.0 & e_pa = 40.0 & e_xcen = 875.0 & e_ycen = 2075.0
       dist_ellipse, mask, imsize, e_xcen, e_ycen, rmajor/rminor, e_pa
       apimage = temporary(apimage) * (mask gt rmajor)
    endif

; ---------------------------------------------------------------------------    
    if (n_elements(file) ne 0L) then if strmatch(file,'*2038_sdss_1_r_c1*') then begin
       rmajor = 730.0 & rminor = 350.0 & e_pa = 35.0 & e_xcen = 760.0 & e_ycen = 1980.0
       dist_ellipse, mask, imsize, e_xcen, e_ycen, rmajor/rminor, e_pa
       apimage = temporary(apimage) * (mask gt rmajor)
    endif

    if (n_elements(file) ne 0L) then if strmatch(file,'*2039_sdss_1_g_c1*') then begin
       rmajor = 730.0 & rminor = 350.0 & e_pa = 35.0 & e_xcen = 760.0 & e_ycen = 1980.0
       dist_ellipse, mask, imsize, e_xcen, e_ycen, rmajor/rminor, e_pa
       apimage = temporary(apimage) * (mask gt rmajor)
    endif

return, apimage
end
    
