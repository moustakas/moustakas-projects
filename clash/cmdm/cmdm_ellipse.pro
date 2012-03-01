;+
; NAME:
; bgt_ellipse
; input:
;   image -- the deblended image containing the galaxy
; optional input:
;   xcen, ycen -- the first guess of the position of the center of the galaxy,
;                 assumed to be in the middle of the image by default
;   na -- the number of grid of semi-major axis
;   amin -- the minimum of the length of the semi-major axis
;           5 pixels by default
; output:
;   ellipse -- a structure containing the fitted parameters
;           xcen
;           ycen
;           sb0
;           ellipticity[na]
;           pa[na]
;           a[na]
;           sb[na]
; bugs:
; COMMENTS:
;    -- the center is constrained within initial guess \pm 5 pixels
; revision history:
;    20-Nov-2008 -- Guangtun Zhu, NYU, started
;    16-Mar-2009 -- Guangtun Zhu, NYU, imivar
;-

pro cmdm_ellipse, image, imivar=imivar, imerr=imerr, ini_xcen0=ini_xcen0, ini_ycen0=ini_ycen0, $
          na=na, amin0=amin0, amax0=amax0, ellipse=ellipse, outimage=outimage, flevel=flevel, $
          nea=nea, limfrac=limfrac, maxiter=maxiter, miniter=miniter, silent=silent, $
          ini_e0=ini_e0, ini_pa0=ini_pa0, status=status, namax=namax

    nx = (size(image))[1]
    ny = (size(image))[2]

    status=1
    print, nx, ny
    if (nx le 4 or ny le 4) then begin
;   if (nx le 10 or ny le 10) then begin
       splog, 'image too small, quiting.....'
       status=-1
       return
    endif

;; xcen0 and ycen0 could be guessed by finding the peaks -- think about it.
;; here assuming it's in the middle of the image
    if (not keyword_set(ini_xcen0)) then ini_xcen0 = float(nx)/2.
    if (not keyword_set(ini_ycen0)) then ini_ycen0 = float(ny)/2.

;;; find the peak of image
;    ii = bsort(image)
;    itop = ceil(0.01*float(nx*ny))
;    iix = ii mod nx
;    iiy = ii/nx
;    xcen_guess = median(float(iix[nx*ny-itop-1L:nx*ny-1L]))
;    ycen_guess = median(float(iiy[nx*ny-itop-1L:nx*ny-1L]))
;    if (not keyword_set(ini_xcen0)) then ini_xcen0 = xcen_guess
;    if (not keyword_set(ini_ycen0)) then ini_ycen0 = ycen_guess

    if (not keyword_set(limfrac)) then limfrac = 4.
    if (not keyword_set(maxiter)) then maxiter = 400L
    if (not keyword_set(miniter)) then miniter = 10L
    if (not keyword_set(namax)) then namax = 300L


    if (not keyword_set(imivar) and keyword_set(imerr)) then imivar = 1./imerr^2
    if (not keyword_set(imivar) and not keyword_set(imerr)) then begin
       splog, 'neither imivar nor imerr is given, quiting ...'
       return
    endif
;;
;;  First guess of parameters. Not working for undeblended version.
;;
    nxcen0 = round(ini_xcen0)
    nycen0 = round(ini_ycen0)
    if nxcen0 lt 4 or nycen0 lt 4 or nxcen0 ge nx-4 or nycen0 ge ny-4 then begin
       status = -1
       return
    endif
    sbcent0 = max(image[nxcen0-4L:nxcen0+4L,nycen0-4L:nycen0+4L])

;; flevel is used to estimate the first guess of pitch angle and semi-axis
    if (not keyword_set(flevel)) then begin
       flevel=[0.05]
    end

;; use the contour at flevel to estimate the first guess
    contour, image, path_xy=xy, path_info=info, closed=0, levels=flevel*max(image), /path_data_coords
;   print, xy[0,*], xy[1,*], total(xy[0,*]*xy[1,*])
;   stop
    if (n_elements(xy) le 4) then begin
       status=-1
       return
    endif
    bgt_sixlin, xy[0,*], xy[1,*], sixa, sixsiga, sixb, sixsigb, status=status
    if (status lt 0) then return
;   if (strmatch(!error_state.MSG, '*SIXLIN*') eq 1) then begin
;      status=-1
;      return
;   endif
    if (not keyword_set(ini_pa0)) then ini_pa0 = atan(sixb[2])

    rr = sqrt((xy[0,*]-ini_xcen0)^2 + (xy[1,*]-ini_ycen0)^2)
    nrr = n_elements(rr)
    rsort = bsort(rr)
    itop = ceil(0.1*nrr)
    rbot = median(rr[rsort[0L:itop-1L]])
    rtop = median(rr[rsort[nrr-itop:nrr-1L]])
    if (not keyword_set(ini_e0)) then ini_e0 = 1. - rbot/rtop

;; Okay, we have the first guess: pa0, e0, xcen0, ycen0

;; 0.396 arcsec/pixel
;; amin
    if (not keyword_set(amin0)) then begin
       amin = 3L ;pixels
    endif else begin
       amin = amin0 > 1L
    endelse
;; da = 1 pixels
;; amax + cen < nx or ny
    
    amax1 = (nxcen0-amin) < (nycen0-amin) < (nx-nxcen0-amin) < (ny-nycen0-amin)
    if (not keyword_set(amax0)) then begin
       amax = amax1
    endif else begin
       amax = amax0 < amax1
    endelse
    amax = long(amax)

;; use sb instead i since i is commonly used for integer
    
    na = amax - amin + 1L
    if na le 5 then begin
      status=-1
       return
    endif

    if na ge namax then begin
       na = namax 
       amax = namax+amin-1L
    endif
    ellipse = {xcen: fltarr(namax), ycen: fltarr(namax), sb0: fltarr(namax), $
        ellipticity: fltarr(namax), pa: fltarr(namax), majora: fltarr(namax), tflux: fltarr(namax), $
        a3: fltarr(namax), b3: fltarr(namax), a4: fltarr(namax), b4: fltarr(namax), $
        a3_err: fltarr(namax), b3_err: fltarr(namax), a4_err: fltarr(namax), b4_err: fltarr(namax), na: na, $
        tflux_ivar: fltarr(namax), sb0_ivar: fltarr(namax), xcenmed: 0.0, ycenmed: 0.0, $
        pafit: fltarr(namax), ellipticityfit: fltarr(namax), $
        a3fit: fltarr(namax), a4fit: fltarr(namax), tfluxfit: fltarr(namax), sb0fit: fltarr(namax), $
        a3fit_err: fltarr(namax), a4fit_err: fltarr(namax), tfluxfit_ivar: fltarr(namax), sb0fit_ivar: fltarr(namax), $
        b3fit: fltarr(namax), b4fit: fltarr(namax), $
        b3fit_err: fltarr(namax), b4fit_err: fltarr(namax), $
        xcenfit: fltarr(namax), ycenfit: fltarr(namax), $
        fit_success:0}

;; First Iteration
;; Are we doing more iterations, e.g., to fix the center?
    for ia = amin, amax do begin
;       print, ia
        xcen0 = ini_xcen0
        ycen0 = ini_ycen0
        e0 = ini_e0
        pa0 = ini_pa0
        a = float(ia)
        niter = 0L
        minmaxp0 = 0.
        ncenbad = 0L
        repeat begin
;; given an image and an ellipse, return the bilinear interpolated flux along the ellipse
;; calculate the derivative of the intensity along the major axis direction
;; simple central finite difference method, f'(x2) = [f(x3)-f(x1)]/[x3-x1]
           sbcoord = bgt_ellipse_sb(image, imivar=imivar, xcen=xcen0, ycen=ycen0, pa=pa0, $
                         ee=e0, majora=a, nea=nea, coords=coords, ea=ea, $
                         sbder=sbder, sbcoord_ivar=sbcoord_ivar)
;; get the fit, ecoeff -- structure
           bgt_ellipse_fit, ea, sbcoord, ecoeff, ini_value=ini_value, fixed=fixed,$
                            sb_ivar=sbcoord_ivar, status=status
           if (status lt 0) then return
           eminus = 1. - e0
           xcen = xcen0 - ecoeff.b1/sbder
           ycen = ycen0 - ecoeff.a1*eminus/sbder
           if (abs(xcen-ini_xcen0) gt 8) then xcen=ini_xcen0
           if (abs(ycen-ini_ycen0) gt 8) then ycen=ini_ycen0
           ee = e0 - 2.*ecoeff.b2*eminus/a/sbder
;          pa = pa0 + 2.*ecoeff.a2*eminus/a/sbder/(eminus^2-1.)
           pa = atan(tan(pa0 + 2.*ecoeff.a2*eminus/a/sbder/(eminus^2-1.)))

           niter = niter + 1L
           
           maxp = abs(ecoeff.a1) > abs(ecoeff.b1) > abs(ecoeff.a2) > abs(ecoeff.b2)
           maxf = maxp/ecoeff.sb0*100.

;; record the initial parameter which generates lowest harmonic amplitude so far
           minmaxp = maxp
           minmaxpf = maxf
           if (niter eq 1L) or (minmaxp le minmaxp0) then begin 
              minmaxp0 = minmaxp
              minmaxpf0 = minmaxpf
              minxcen = xcen0
              minycen = ycen0
              minee = e0
              minpa = pa0
              minsb = ecoeff.sb0
           endif

           if (ee le 0. or ee gt 1.0) then ee = ini_e0
           if (xcen-a-1. le 0. or ycen-a-1. le 0. or xcen+a+1. ge  nx or ycen+a+1. ge ny)  then begin
              if (not keyword_set(silent)) then begin
                 print, ini_xcen0, xcen, nx
                 print, ini_ycen0, ycen, ny
              endif
              xcen = ini_xcen0
              ycen = ini_ycen0
           endif
           xcen0 = xcen
           ycen0 = ycen
           e0 = ee
           pa0 = pa
           sb0 = ecoeff.sb0

        endrep until ((niter ge miniter) and ((niter eq maxiter) or (maxf le limfrac)))
        
        if ((niter le maxiter) or (maxf le limfrac)) then begin
           xcen0 = minxcen
           ycen0 = minycen
           e0 = minee
           pa0 = minpa
           sb0 = minsb
        endif
        
        ii = ia - amin
        ellipse.xcen[ii] = xcen0
        ellipse.ycen[ii] = ycen0
        ellipse.sb0[ii] = sb0
        ellipse.ellipticity[ii] = e0
        ellipse.pa[ii] = pa0
        ellipse.majora[ii] = a

;; total flux within the ellipse
        sbcoord = bgt_ellipse_sb(image, imivar=imivar, xcen=xcen0, ycen=ycen0, pa=pa0, ee=e0, $
                              majora=a, nea=nea, coords=coords, ea=ea, $
                              sbder=sbder, sbcoord_ivar=sbcoord_ivar)

        ellipse.sb0_ivar[ii] = median(sbcoord_ivar)
        sub_1d = polyfillv(coords[0,*], coords[1,*], nx, ny)
        if (n_elements(sub_1d) gt 1L) then begin
           sub_1dii = where(finite(imivar[sub_1d]) ne 0, mm)
           if (mm gt 0) then begin 
              ellipse.tflux[ii] = total(image[sub_1d[sub_1dii]])
              ellipse.tflux_ivar[ii] = 1./(total(1./imivar[sub_1d[sub_1dii]]))
           endif
        endif
;; fit a3, a4, b3, b4
        bgt_ellipse_fit2, ea, sbcoord, ecoeff2, $
                ini_value=ini_value, fixed=fixed, $
                sb_ivar=sbcoord_ivar
        ellipse.a3[ii] = ecoeff2.a3
        ellipse.a3_err[ii] = ecoeff2.a3_err
        ellipse.b3[ii] = ecoeff2.b3
        ellipse.b3_err[ii] = ecoeff2.b3_err
        ellipse.a4[ii] = ecoeff2.a4
        ellipse.a4_err[ii] = ecoeff2.a4_err
        ellipse.b4[ii] = ecoeff2.b4
        ellipse.b4_err[ii] = ecoeff2.b4_err
        if (not keyword_set(silent)) then begin
           print, xcen0, ycen0, ini_xcen0, ini_ycen0
        endif
        
    endfor

;; to-do list
;; fit pa, ee with a low order polynomial, return, xcenmedian, ycenmedian
;; recalculate tflux and sb0

;  if (na lt 5) then fit_success=0
   if (na ge 5) then begin
      ellipse.fit_success=1
      npoly = 3
      nsigma = 3.

;; don't use the outer 10% to fit the profile
      nnaa  = round(na*0.90)
      aa = fltarr(npoly+1L,na)
      for i = 0L, npoly do begin
         aa[i,*]=findgen(na)^i
      endfor

      ivar = fltarr(na)+1.
      infnan = where(finite(ellipse.sb0_ivar[0:na-1]) eq 0, nin)
      if (nin gt 0L) then ivar[infnan] = 0.
      if (nnaa lt na) then ivar[nnaa-1:na-1] = 0.

      yy = ellipse.pa[0:na-1]
      hogg_iter_linfit, aa, yy, ivar, coeffs, nsigma=nsigma
;     ellipse.pafit[0:na-1] = aa##coeffs
      ellipse.pafit[0:na-1] = smooth(ellipse.pa[0:na-1],5)

      yy = ellipse.ellipticity[0:na-1]
      hogg_iter_linfit, aa, yy, ivar, coeffs, nsigma=nsigma
;     ellipse.ellipticityfit[0:na-1] = aa##coeffs
      ellipse.ellipticityfit[0:na-1] = smooth(ellipse.ellipticity[0:na-1],5)

      xcenall = ellipse.xcen[0:na-1]
;     yy = ellipse.xcen[0:na-1]
;     hogg_iter_linfit, aa, yy, ivar, coeffs, nsigma=nsigma
;     ellipse.xcenfit[0:na-1] = aa##coeffs
      ellipse.xcenfit[0:na-1] = smooth(ellipse.xcen[0:na-1],5)

      ycenall = ellipse.ycen[0:na-1]
;     yy = ellipse.ycen[0:na-1]
;     hogg_iter_linfit, aa, yy, ivar, coeffs, nsigma=nsigma
;     ellipse.ycenfit[0:na-1] = aa##coeffs
      ellipse.ycenfit[0:na-1] = smooth(ellipse.ycen[0:na-1],5)

      ii = where(ivar gt 0., nn)
      if (nn gt 1) then begin 
         ellipse.xcenmed = median(xcenall[ii])
         ellipse.ycenmed = median(ycenall[ii])
      endif

      bgt_ellipse_allfit, image, ellipse, imivar=imivar, amin=amin, amax=amax
;; total flux within the ellipse
;     for ia = amin, amax do begin
;         ii = ia - amin
;; float (ia) should be the same as ellipse.majora[ii]
;         a = float(ia)
;         xcen0 = ellipse.xcenfit[ii]
;         ycen0 = ellipse.ycenfit[ii]
;         pa0 = ellipse.pafit[ii]
;         e0 = ellipse.ellipticityfit[ii]
;         sbcoord = bgt_ellipse_sb(image, imivar=imivar, xcen=xcen0, ycen=ycen0, pa=pa0, ee=e0, $
;                   majora=a, nea=nea, coords=coords, ea=ea, $
;                   sbder=sbder, sbcoord_ivar=sbcoord_ivar)

;         ellipse.sb0fit_ivar[ii] = median(sbcoord_ivar)
;         sub_1d = polyfillv(coords[0,*], coords[1,*], nx, ny)
;         if (n_elements(sub_1d) gt 1L) then begin
;            sub_1dii = where(finite(imivar[sub_1d]) ne 0, mm)
;            if (mm gt 0) then begin
;               ellipse.tfluxfit[ii] = total(image[sub_1d[sub_1dii]])
;               ellipse.tfluxfit_ivar[ii] = 1./(total(1./imivar[sub_1d[sub_1dii]]))
;            endif
;         endif
;; fit a3, a4, b3, b4
;         bgt_ellipse_fit2, ea, sbcoord, ecoeff2, $
;                 ini_value=ini_value, fixed=fixed, $
;                 sb_ivar=sbcoord_ivar
;         ellipse.a3fit[ii] = ecoeff2.a3
;         ellipse.a3fit_err[ii] = ecoeff2.a3_err
;         ellipse.b3fit[ii] = ecoeff2.b3
;         ellipse.b3fit_err[ii] = ecoeff2.b3_err
;         ellipse.a4fit[ii] = ecoeff2.a4
;         ellipse.a4fit_err[ii] = ecoeff2.a4_err
;         ellipse.b4fit[ii] = ecoeff2.b4
;         ellipse.b4fit_err[ii] = ecoeff2.b4_err
;         ellipse.sb0fit[ii] = ecoeff2.sb0
;         if (not keyword_set(silent)) then begin
;            print, xcen0, ycen0, ini_xcen0, ini_ycen0
;         endif
;     endfor

   endif
   
   infnan = where(finite(ellipse.sb0_ivar[0:na-1]) eq 0, nin)
   if (nin gt 0L) then ellipse.sb0_ivar[infnan] = 0.

end
