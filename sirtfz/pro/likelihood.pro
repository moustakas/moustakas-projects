pro likelihood, oflux, oflux_error, mflux, chi2_matrix, lhood_matrix, constant_matrix
;+
; NAME:
;	LIKELIHOOD
;
; PURPOSE:
;	Evaluate the likelihood surface that compares the data and the
;	model fluxes for photometric redshifts.
;
; CALLING SEQUENCE:
;	likelihood, oflux, oflux_error, mflux, chi2_matrix, $
;	lhood_matrix, constant_matrix
;
; INPUTS:
;	oflux       - observed flux array [nbands]
;	oflux_error - observed flux error array [nbands]
;	mflux       - model fluxes [nseds,nz,nbands]
;
; OUTPUTS:
;	chi2_matrix     - chi-squared surface [nseds,nz]
;	lhood_matrix    - likelihood surface [nseds,nz]
;	constant_matrix - ratio of observed to model fluxes surface
;                         [nseds,nz]
;
; COMMENTS:
;	See Benitez 2000 for the analytic expressions for the chi2. 
;
; MODIFICATION HISTORY:
;	J. Moustakas, 2001 Jan 15, U of A, written
;	jm01jul30uofa, documented
;	jm01sep12uofa, minor modifications
;-

; create the arrays we will need

    msize = size(mflux,/dimension)
    nseds = msize[0]
    nz = msize[1]

    f_tt = dblarr(nseds,nz)
    f_0t = dblarr(nseds,nz)
    
; fill the arrays (benitez 2000, equation 9)
    
    f_00 = total((oflux/oflux_error)^2.0D)
    for i = 0L, nseds-1L do begin
       f_tt[i,*] = (reform(mflux[i,*,*]))^2.0D # (1.0D/oflux_error^2.0D)
       f_0t[i,*] = reform(mflux[i,*,*]) # (oflux/oflux_error^2.0D)
    endfor

; evaluate chi-squared and the likelihood (benitez 2000, equation 7)

    chi2_matrix = f_00 - f_0t^2.0D / f_tt
    bad = where(finite(chi2_matrix) eq 0B,nbad)

    lhood_matrix = exp(-0.5D*chi2_matrix)
    if nbad ne 0L then lhood_matrix[bad] = 0.0 ; set infinities to zero
    norm = total(lhood_matrix) 
    if norm gt float(0) then lhood_matrix = lhood_matrix/norm ; normalize

    constant_matrix = f_0t / f_tt ; constant that minimizes chi2
    
return
end
    
