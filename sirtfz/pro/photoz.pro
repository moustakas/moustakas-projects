pro photoz, oflux, oflux_error, mflux, zarray, pdz=pdz, pzstruc
;+
; NAME:
;	PHOTOZ
;
; PURPOSE:
;	Calculate the photometric redshift of a galaxy.
;
; CALLING SEQUENCE:
;	photoz, oflux, oflux_error, mflux, zarray, [pdz=], pzstruc
;
; INPUTS:
;	oflux       - observed flux [nbands] 
;	oflux_error - flux error [nbands] 
;	mflux       - model flux cube [nseds,nz,nbands]
;	zarray      - redshift array [nz]
;
; OPTIONAL INPUTS:
;	pdz          - confidence interval [0,1] within which to compute
;                      the redshift uncertainty [default 0.95]
; OUTPUTS:
;	pzstruc     - structure detailing the parameters of the fit
;
; COMMENTS:
;
; PROCEDURES USED:
;	LIKELIHOOD, CUMULATE
;
; MODIFICATION HISTORY:
;	J. Moustakas, 2001 May 13, U of A
;	jm01aug2uofa, documented
;	jm01sep14uofa, added mean and median photo-z estimation and an
;	estimate on the uncertainty based on the cumulative PDF
;-
    
    if keyword_set(pdz) then pdz = (pdz > 0.0) < 1.0 else pdz = 0.95

    msize = size(mflux,/dimension)
    nseds = msize[0]
    nz = msize[1]
    nbands = msize[2]

; evaluate chi2

    likelihood, oflux, oflux_error, mflux, chi2_matrix, lhood_matrix, constant_matrix

    maxprob = max(lhood_matrix,mindx)
    if maxprob eq float(0) then minchi2 = min(chi2_matrix,mindx)

; ----------------------------------------------------------------------
; compute the mode of the distribution
; ----------------------------------------------------------------------
    
    zindx = mindx/nseds     ; redshift index
    tindx = mindx mod nseds ; SED type index

    z_mode = zarray[zindx]

; ----------------------------------------------------------------------
; compute the mean (probability-weighted centroid of the likelihood
; surface)
; ----------------------------------------------------------------------

    zpdf = total(lhood_matrix,1) ; sum over the type dimension
    zarea = total(zpdf)
    if zarea gt float(0) then z_mean = total(zpdf*zarray)/zarea else z_mean = 0.0
;   zindx_mean = interpol(findgen(nz),zarray,z_mean)
    
; ----------------------------------------------------------------------
; compute the median using the cumulative distribution function
; ----------------------------------------------------------------------

    if zarea gt float(0) then begin

       z_cumulative = cumulate(zpdf)
       z_median = interpol(zarray,z_cumulative,0.5)
;      zindx_median = interpol(findgen(nz),zarray,z_median)

       if z_median lt float(0) then begin

          z_median = 0.0
          sigmaz = [0.0,0.0]

       endif else begin
       
; pdz percent of p(z) is contained in the interval [zmin,zmax]
    
          zmin = interpol(zarray,z_cumulative,0.5-pdz/2.0)
          zmax = interpol(zarray,z_cumulative,0.5+pdz/2.0)
          sigmaz = [z_median-zmin,zmax-z_median]

       endelse

    endif else begin

       z_median = 0.0
       sigmaz = [0.0,0.0]

    endelse
       
; ----------------------------------------------------------------------
; find the minimum of chi2 (using the mode)

;   chi2 = interpolate(chi2_matrix,tindx,zindx_mean)
    
    chi2 = chi2_matrix[tindx,zindx]    ; chi2
    chi2_nu = chi2 / float(nbands-1.0) ; chi2_nu (reduced chi2)
    
    constant = constant_matrix[tindx,zindx] ; amplitude
    
    pzstruc.zphot_mode = float(z_mode)
    pzstruc.zphot_mean = float(z_mean)
    pzstruc.zphot_median = float(z_median)
    pzstruc.sigmaz = float(sigmaz)
    pzstruc.chi2_nu = float(chi2_nu)
    pzstruc.zpdf = zpdf
    pzstruc.constant = float(constant)
    pzstruc.zindx = long(zindx)
    pzstruc.tindx = long(tindx)
    
return
end
