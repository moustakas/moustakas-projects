;**************************
pro getwave,h,wave
;+
; PURPOSE:
;	Create a wavelength vector from a FITS header of a spectrum
;	The FITS header must contain a CRVAL1 keyword and a CRPIX1 or
;	CD1_1 keyword.
;
; CALLING SEQUENCE:
;	GETWAVE, hdr, wave
;
; INPUTS:
;	hdr - FITS header for a 1-d spectrum, string array
;
; OUTPUTS:
;	wave - Wavelength vector.    
;-

 zparcheck, 'GETWAVE', h, 1, 7, 1, 'FITS header'

 xsize = sxpar(h,'NAXIS1')
 if xsize LE 0 then message, $
	'ERROR - FITS header does not contain positive value of NAXIS1'
 w0 = sxpar(h,'CRVAL1')
 if !ERR EQ -1 then message,'FITS header does not contain CRVAL1 value'
 wdelt = sxpar(h,'CDELT1')
 if !ERR EQ -1 then begin
	wdelt = sxpar(h,'CD1_1')
        if !ERR EQ -1 then message, $
             'FITS header does not contain CDELT1 or CD1_1 value'
 endif
 
 wave = w0 + wdelt*findgen(xsize)

 return
 end

;+
; NAME:
;      remove_telluric
;
; PURPOSE:
;      remove effects of telluric absorption bands in the 1D spectra
;
; CALLING SEQUENCE:
;      remove_telluric, spec1d, airmass, plot=plot
; 
; INPUTS:
;      spec1d - ss1d structure containing spec, lambda, ivar
;      airmass - airmass of spec1d object
;
; KEYWORD PARAMETERS:
;      plot - if set, splot output of before and after spectrum
;
; OUTPUTS:
;      returns spec1d structure with corrected .spec spectrum
;
;
; MODIFICATION HISTORY:
;   alc 03oct02
;   md  070ct02 - use common block to store processed bands
;   mcc 05nov02 - minor changes to make routine more robust
;   md  07feb03 - deal with grating variations requiring different smoothing
;   ogt 16may12 - enable handling of 600 line/mm data
;-

;**************************
pro read_aband,aa
 
  aband_dir = getenv('DEEP2_SPEC2D')+'/etc/'
;  aband_dir = getenv('IDLSPEC2D_DIR')+'/etc/'
  a_flux=readfits(aband_dir +'aband.fits',heada, /silent)
  getwave,heada,a_lambda
 ; normalize spectrum
  a_fluxn=float(a_flux)/((total(a_flux[1400:1499])+total(a_flux[2500:2599]))/200.)



  minab=1510
  maxab=2110

; deweight a band itself in continuum fit
  errors=a_lambda*0+1.
  errors[minab:maxab]=1.E10

 ; remove continuum shape - quadratic polynomial fit -- LINEAR WAS BETTER!!!
; jm14jun02siena - use the much stabler linfit()
   aacoeff = linfit(a_lambda, a_fluxn,measure = errors, yfit=aband_fit) 
;  aacoeff = svdfit(a_lambda, a_fluxn, measure=errors,2, yfit=aband_fit) 
   ; shift fit up to match continuum of a_flux - NO LONGER REQUIRED
;   shift = mean(a_fluxn[1000:1200])-mean(aband_fit[1000:1200])
;print,shift
;   aband_fit = aband_fit+shift
   a_flux = a_fluxn/ aband_fit 

   aa = {a_lambda:float(a_lambda[minab:maxab]), $
         a_flux:float(a_flux[minab:maxab]) }
   
return
end
;**************************
pro read_bband, bb
   bband_dir = getenv('DEEP2_SPEC2D')+'/etc/'
;  bband_dir = getenv('IDLSPEC2D_DIR')+'/etc/'
   b_flux=readfits(bband_dir +'bband.fits',headb, /silent)
   getwave,headb,b_lambda
; avoid the false feature (absorption spike) at 6556 A.

   dex1 = min(where(b_lambda gt 6530.))
   dex2 = min(where(b_lambda gt 6580.))

; don't fit continuum around it or the B band itself:

   minbb = min(where(b_lambda gt 6840.))
   maxbb = min(where(b_lambda gt 6960.))
   minuseful = min(where(b_lambda gt 6650.))

   errors=b_lambda*0+1.
   errors[minbb:maxbb]=1.E10
   errors[dex1:dex2]=1.E10
   errors[0:minuseful] = 1.E20

 ; normalize spectrum
   b_flux=float(b_flux)/((total(b_flux[minbb-99:minbb])+ $
                          total(b_flux[maxbb:maxbb+99])) /200.)

; trim spectrum - BETTER TO DEWEIGHT BAD REGIONS AS ABOVE
;   b_lambda = b_lambda[300:3200] & b_flux = b_flux[300:3200]

 ; remove continuum shape - quadratic polynomial fit - LINEAR WAS BETTER
; jm14jun02siena - use the much stabler linfit()
   bbcoeff = linfit(b_lambda, b_flux,measure = errors, yfit=bband_fit) 
;  bbcoeff = svdfit(b_lambda, b_flux,measure = errors, 2, yfit=bband_fit) 

   ; shift fit up to match continuum of b_flux - NOT NECESSARY!
;   shift = mean(b_flux[1000:1200])-mean(bband_fit[1000:1200])
;   bband_fit = bband_fit+shift
   b_flux = b_flux/ bband_fit 

; I HAVE NO IDEA WHY YOU USED SUCH A LARGE WINDOW BEFORE.  THE ONLY
;   SIGNIFICANT FEATURES IN THIS WINDOW IN CHUCK'S TEMPLATE ARE _NOT_
;   IN GEOFF MARCY'S.

   bb = {b_lambda:float(b_lambda[minbb:maxbb]), $
         b_flux:float(b_flux[minbb:maxbb]) }

return
end
;**************************
; used to get a pixel shift between wavelength array of the input
;     spectrum compared to a template spectrum

pro get_shift, in_lambda, in_flux, temp_lambda, temp_flux,  pixel_shift

;apodize the ends!
   npix = n_elements(in_flux) 
   for i=0, npix/10-1 do begin 
       in_flux[i] = in_flux[i]*sin(i*!pi/2/(npix/10))
       in_flux[npix-1-i] = in_flux[npix-1-i]*sin(i*!pi/2/(npix/10))
   endfor

   npix = n_elements(temp_flux) 
   for i=0, npix/10-1 do begin 
       temp_flux[i] = temp_flux[i]*sin(i*!pi/2/(npix/10))
       temp_flux[npix-1-i] = temp_flux[npix-1-i]*sin(i*!pi/2/(npix/10))
   endfor

; restrict to +/- 7 pixel shift
   ncor = 15 
   noff = fix(ncor/2)
   xcor=fltarr(ncor)
   xx=findgen(ncor)-noff

   npix = n_elements(in_lambda) < n_elements(temp_lambda) ; shorter array
   xnorm = sqrt(variance(temp_flux)*variance(in_flux))    ; normalization
; do the cross-correlation
   for j=0,ncor-1 do xcor[j]=total(temp_flux*shift(in_flux,j-noff))/npix/xnorm
   xpeak = max(xcor, ipeak) ;peak of x-corr.
   xxs = xx[ipeak-2 > 0:ipeak+2 < ncor-1]
   xcors = xcor[ipeak-2 > 0:ipeak+2 < ncor-1] ;select subset around peak
   fit=poly_fit(xxs,xcors,2)
   pixel_shift=-.5*fit[1]/fit[2] ; get central value of fit

;window,1
;plot,xcor
return
end

;**************************


pro remove_telluric, spec1d, airmass, plot = plot


   common comm_telluric, aaa,  aspl, aband_airmass, $
                    bbb,  bspl, bband_airmass, smoothed, bsmoothed ;band data

   flux = spec1d.spec
   lambda = spec1d.lambda
   ivar = spec1d.ivar
   npixsp = n_elements(lambda)

   if n_elements(smoothed) eq 0 then smoothed = 0
   smoothed=0
   bsmoothed = 1

; define the angstrom shift.
   angstrom_shift = -0.77

   if keyword_set(aaa) eq 0 then begin ;input band information, if needed
;;;; DO A-BAND FIRST
;get aband template and wavelength
     read_aband, aa
     aaa = aa
; correct template for airmass of observed object
     aband_airmass = 1.0
   endif else aa = aaa




   npix = n_elements(aaa.a_flux)

;  get ratio of wavelength scales for possible smoothing of feature
   sigma = (lambda[1001]-lambda[1000])/(aa.a_lambda[1]-aa.a_lambda[0])
   smooth = (sigma gt 1.1) ? 1:0
   if smooth then begin ;gaussian smooth by ratio of grating grooves
     kernel = fltarr(15)
     for i=0, 14 do kernel[i] = exp(-.5*(i-7.)^2/sigma^2)
     kernel = kernel/total(kernel)
     aa.a_flux = convol(aa.a_flux, kernel, /center, /edge_truncate)
 endif



; generate spline interpolation to match data spectral dispersion 
   aspl = spl_init(aa.a_lambda, aa.a_flux)

; trim input spectrum to just cover A band
   a_lambda_low = where(lambda lt aa.a_lambda[0], lowcnt)
   if lowcnt eq 0 then almin = 0 $
      else almin=max(a_lambda_low)+1
   npix = n_elements(aa.a_flux)
   a_lambdamax = where(lambda gt aa.a_lambda[npix-1], upcnt)
   if upcnt eq 0 then almax = npixsp-1 $
      else almax = min(a_lambdamax)-1
; check that A band falls on the 1-d spectrum...
   inrange = where(lambda gt aa.a_lambda[0] and $
                   lambda lt aa.a_lambda[npix-1], acnt)
   if acnt eq 0 then begin
       print, 'A-band not in spectral window...' 
   endif else begin 
       if almax lt almin then print, '(remove_telluric.pro) error with wavelength solution!' $
       else begin

; interpolate A band template to same wavelength array as input
; spectrum and scale according to airmass
           aband_flux = spl_interp(aa.a_lambda, aa.a_flux, aspl, $
                          lambda[almin:almax])^((airmass/aband_airmass)^0.55)


           pixel_shift=0.




   ; now shift the template to match the input spectrum
           aa.a_lambda = aa.a_lambda-angstrom_shift

   ; re-trim the input spectrum to just cover the shifted A band
           a_lambda_low = where(lambda lt aa.a_lambda[0], lowcnt) 
           if lowcnt eq 0 then almin = 0 $
           else almin = max(a_lambda_low)+1
           a_lambdamax = where(lambda gt aa.a_lambda[npix-1], upcnt)
           if upcnt eq 0 then almax = npixsp-1 $
           else almax = min(a_lambdamax)-1
       
    ; now re-interpolate with the new wavelength scale
           aspl = spl_init(aa.a_lambda, aa.a_flux)


    ; new A band flux
           aband_flux = spl_interp(aa.a_lambda, aa.a_flux, aspl, $
                          lambda[almin:almax])^((airmass/aband_airmass)^0.55)

           aband_lambda = lambda[almin:almax]


;FOR LOW-RESOLUTION DATA 
           ; test whether low- or high-resolution
           deriv=mean(deriv(spec1d.lambda))
           if deriv gt 0.5 and smoothed eq 0 then begin
               ; scale template and construct kernel for convolution
               aband_flux=aband_flux^1.15 
               x=(dindgen(21) - 10.)*0.66
               kernel2=exp(-(x^2)/(2*1.19^2))
               kernel_norm=kernel2/int_tabulated(x,kernel2)
               
               ; new aband_flux for low-res data
               new_aband_flux=convol(aband_flux, kernel_norm, /CENTER)
               aband_flux=new_aband_flux/max(new_aband_flux)     ;normalize to 1
               aband_flux[where(aband_flux eq 0.)] = 1.      ;set ends to 1

		smoothed = 1
		bsmoothed = 0
            endif


    ; correct flux and ivar

           ;splot,lambda,spec1d.spec,psym=10

           spec1d.spec[almin:almax] = spec1d.spec[almin:almax]/aband_flux 
           spec1d.ivar[almin:almax] = spec1d.ivar[almin:almax]*aband_flux^2 
           
       endelse
   endelse
    
;;; NOW DO B-BAND
;get aband template and wavelength
   if keyword_set(bbb) eq 0 then begin ;input band information, if needed
       read_bband, bb
       bbb = bb
       bband_airmass = 1.0
   endif else bb = bbb


   npix = n_elements(bb.b_flux)

  if smooth then $ ;gaussian smooth by grating ratio, if needed
     bb.b_flux = convol(bb.b_flux, kernel, /center, /edge_truncate)


; shift the template to match the input spectrum - shift by same as A band
   bb.b_lambda = bb.b_lambda-angstrom_shift

; generate spline interpolation to match data spectral dispersion 
   bspl = spl_init(bb.b_lambda, bb.b_flux)

; trim input spectrum to just cover B band
   b_lambda_low = where(lambda lt bb.b_lambda[0], lowcnt) 
   if lowcnt eq 0 then blmin = 0 else blmin = max(b_lambda_low)
   npix = n_elements(bb.b_lambda)
   b_lambdamax = where(lambda gt bb.b_lambda[npix-1], upcnt)
   if upcnt eq 0 then blmax = npixsp-1 else blmax = min(b_lambdamax)

; check that B band falls on the 1-d spectrum...
   inrange = where(lambda gt bb.b_lambda[0] and $
                   lambda lt bb.b_lambda[npix-1], bcnt)
   if bcnt eq 0 then print, 'B-band not in spectral window...' $
   else begin 
       if blmax lt blmin then print, '(remove_telluric.pro) error with wavelength solution!' $
       else begin

; interpolate B band template to same wavelength array as input
           bband_flux = spl_interp(bb.b_lambda, bb.b_flux, bspl, $
                                   lambda[blmin:blmax])^((airmass/bband_airmass)^0.55)

;FOR LOW-RESOLUTION DATA 
           ; test whether low- or high-resolution
           ;deriv=mean(deriv(spec1d.lambda))
           if deriv gt 0.5 and bsmoothed eq 0 then begin
               ; scale template and construct kernel for convolution
               bband_flux=bband_flux^1.15

               ; the following lines are repeated from A-band section above for reference 
               ;x=(dindgen(21) - 10.)*0.66 
               ;kernel2=exp(-(x^2)/(2*1.19^2))
               ;kernel_norm=kernel2/int_tabulated(x,kernel2)
               
               ; new aband_flux for low-res data
               new_bband_flux=convol(bband_flux, kernel_norm, /CENTER)
               bband_flux=new_bband_flux/max(new_bband_flux)     ;normalize to 1
               bband_flux[where(bband_flux eq 0.)] = 1.     ;set ends to 1

            endif

           
; correct flux and ivar
           spec1d.spec[blmin:blmax] = spec1d.spec[blmin:blmax]/bband_flux 
           spec1d.ivar[blmin:blmax] = spec1d.ivar[blmin:blmax]*bband_flux^2
           ;soplot,lambda,spec1d.spec,color=1,psym=10
       endelse
   endelse

   if keyword_set(plot) then begin
       window,0
       colortable1
       plot, lambda, flux, xr = [7550, 7720],color=1
       oplot, lambda, spec1d.spec, color = 2
       wait, 5
       window,1
       plot,lambda,flux,color=1,xr=[6700,7000]
       oplot,lambda,spec1d.spec,color=2
   endif


   return
end





