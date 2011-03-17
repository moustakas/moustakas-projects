pro boot_response, iter, lf, bootmags, booterrs, lfstart, lfend, trgbmag, log=log
;+
; NAME: 
;	BOOT_RESPONSE
;
; PURPOSE:
;	Bootstrap resample a luminosity function, the response
;	(derivative) in each magnitude bin, and the sigma amplitude of
;	the response in each bin.
;
; INPUTS:
;	lfstart, lfend : search for the TRGB between these two *indx*
;	elements 
;
; KEYWORD PARAMETERS:
;	log :	bootstrap resample a logarithmic luminosity function
;	smf :	generate a Sakai, Madore & Freedman luminosity
;		function; set this keyword equal to the desired
;		"smooth factor" 
;
; OUTPUTS:
;	respsig :	standard deviation of the response in each
;			luminosity function magnitude bin. 
;	trgbmag :	magnitude of the maximum response (TRGB)
;
; COMMON BLOCKS:
;
; PROCEDURE:
;
; MODIFICATION HISTORY:
;	John Moustakas, 2000 June 26, UCB
;-

	boothist = fltarr(lf.nhist,iter) ; array of resampled LF's
        trgbmag = fltarr(iter)  	 ; magnitude of the TRGB (maximum response)

        for i = 0L, iter-1L do boothist[*,i] = $
          histogram(bootmags[*,i],bin=lf.binsize,min=lf.minmag,max=lf.maxmag)

; calculate the response
        
        if keyword_set(log) then $
          bootresp = shift(alog10(float(boothist)>1.),-1,0)-shift(alog10(float(boothist)>1.),1,0) else $
          bootresp = shift(boothist,-1,0)-shift(boothist,1,0)

; constrain the endpoints

        bootresp[0L:1L] = 0. & bootresp[lf.nhist-2L:lf.nhist-1L] = 0.

        for j = 0L, iter-1L do begin ; find the TRGB from the maximum response

            linsboot = smooth2(float(boothist[*,j]),4) ; Gaussian smooth the linear histogram
            bootnoise = sqrt(2.*(linsboot>1.)) 	  ; Poisson noise per bin
            
            if keyword_set(log) then $
              maxresp = max(bootresp[lfstart:lfend,j]*bootnoise[lfstart:lfend],trgb) else $ ; log
              maxresp = max(bootresp[lfstart:lfend,j]/bootnoise[lfstart:lfend],trgb) ; linear

            trgbmag[j] = (lf.mag[lfstart:lfend])[trgb]

        endfor

return
end
