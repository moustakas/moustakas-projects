FUNCTION readfilt, filter, vega=vega, AB=AB
; kf=readfilt('bessell_k_001.dat')   
    filtdir =  '/data/orfasay/leonidas/Templates/Filters/nonhst/'
   
;   readcol,filtdir+'johnson_b_001.dat',bjohnsl,bjohnst
;   readcol,filtdir+'johnson_v_001.dat',vjohnsl,vjohnst
;   readcol,filtdir+'cousins_r_001.dat',rcousl,rcoust
;   readcol,filtdir+'cousins_i_001.dat',icousl,icoust
;   readcol,filtdir+'bessell_j_001.dat',jbessl,jbesst
;   readcol,filtdir+'bessell_h_001.dat',hbessl,hbesst
;   readcol,filtdir+'bessell_k_001.dat',kbessl,kbesst
   
   IF n_elements(vega) NE 0 THEN BEGIN 
      A0dir = '/home/ioannis/sirtf/data/'
      readcol,A0dir+'aov_kur_bb.sed',a0w,a0f
      return,transpose([[a0w],[a0f]])
   ENDIF 

   readcol,filtdir+filter,filtw,filtf
   return,transpose([[filtw],[filtf]])
   
END 
;--
FUNCTION filtflux,sedw,sedf,farr,z=z
   
   IF n_elements(z) EQ 0 THEN z = 0.

; match the sed wavelength scale to the filter's, and integrate
;   intsedfilt = farr[1,*] * interpol(sedf,sedw,farr[0,*]/(1.+z))
;   intsedflux = int_tabulated(farr[0,*],intsedfilt) / (1.+z)

; match the filter wavelength scale to the sed's, and integrate   
   widx = where(sedw GT min(farr[0,*]) AND sedw LT max(farr[0,*]))
   intsedfilt = interpol(farr[1,*],farr[0,*]/(1.+z),sedw[widx])
   intsedflux = int_tabulated(sedw[widx],intsedfilt*sedf[widx]) / (1.+z)

   return,intsedflux
   
END
;--
FUNCTION zeromag,farr
;   return zeropoint for m = 0 wrt A0V

; Scale total flux (in Lo =  solar) to obtain vmag=0
; vmag = 2.422-2.5*alog10(fv)
   
   a0 = readfilt(/vega)
   vf = readfilt('johnson_v_001.dat')
   f1 = filtflux(a0[0,*],a0[1,*],vf)
   f2 = 10.d0^(0.4*2.422)
   rr = f2/f1
   a0[1,*] = rr*a0[1,*]
   zeromag = 2.5*alog10(filtflux(a0[0,*],a0[1,*],farr))
   return,zeromag
   
END 
;--
FUNCTION filtcolor,sedw,sedf,farr1,farr2,z=z
   return,filtflux(sedw,sedf,farr1,z=z)-filtflux(sedw,sedf,farr2,z=z)
END 
;--
FUNCTION filterflux,filter,sedlam,sed,z
; To convert f_lambda spectrum to f_nu:  fnu = flam * lam * lam / 2.99793e18
; To convert f_nu spectrum to f_lambda:  flam = 2.99793e18 * f_nu / lam / lam
; To convert f_nu spectrum to AB magnitudes: mab = -2.5 * alog10(fnu) - 48.594
; To convert AB magnitudes to f_nu spectrum: fnu = 10^(-0.4 * (mab+48.594))
; set Uo_flux = 1320   # Jky  
; set Bo_flux = 4380   # Jky  
; set Vo_flux = 3490   # Jky  
; set Io_flux = 2310   # Jky  
; set Jo_flux = 1600   # Jky  
; set Ho_flux = 1080   # Jky  
; set Ko_flux = 680    # Jky  
; Compute V magnitude for a 1 Mo galaxy
; It is -27.5 magnitudes brighter for a 1E11 Mo galaxy
; vmag=2.422-2.5*alog10(fx(nv))
; abset 1 # Convert GISSEL flambda SED to ab mags, normalized at m(5556)=20.0
   sol = 2.99792458d10          ; cm/s
END

