function qe_correction2,spec1d, header, params, paramsendr, lambda,throughput_table,ncorrect=ncorrect,telluric=telluric
;+
; NAME:
;      QE_CORRECTION2
;
;
; PURPOSE:
;      Corrects the counts per hour (SPEC) in a spec1d structure for chip-to-chip
;      variations and throughput. QE_CORRECTION2
;      uses pre-tabulated parameters for the chip-to-chip variation
;      and red-end throughput derived from coadded spectra, augmenting Ricardo Schiavon's
;      throughput measurements.  Loosely based on a code from Renbin
;      Yan, but using new measurements.
;       
; SYNTAX
;      spec1d_corr=qe_correction2(spec1d,header,[params, paramsendr, lambda,$
;               throughput_table, /TELLURIC,NCORRECT=ncorrect])
;
;
; CATEGORY:
;      spec1d 
;             
;
; INPUTS:
;      spec1d = a 1d spectrum structure, as read in by FILL_GAP. 
;      header = the FITS header corresponding to spec1d; can be read in via
;        FILL_GAP at the same time as spec1d using the HEADER keyword
;
;      OUTPUTS: 
; 
;      spec1d_corr=version of the spec1d structure, with SPEC and IVAR
;         corrected for QE variations and the
;         instrumental/telescope/atmosphere response curve from
;         Ricardo (modified at the red end)
;
; OPTIONAL INPUTS:
;
;      params = the table of parameters used to correct the quantum
;         efficiencies of chips 1,2,4,5,6,8 to match chips 3/7, as well
;         as a small correction for the chip3/7 jump.  Automatically read
;         in if not provided or the routine is being run for the first time.
;
;      paramsendr = the parameters found to correct for the fall off
;          of throughput at the red end of the
;          spectra; Ricardo's correction appears to be
;          insufficient.  Automatically read
;          in if not provided or the routine is being run for
;          the first time.
;      lambda   = wavelength array for Ricardo's throughput table
;      throughput_table = throughput array for Ricardo's
;          throughput table.  Automatically read
;          in if not provided or the routine is being run 
;          for the first time.
;
;  KEYWORD INPUTS
;      /TELLURIC = set this keyword to have QE_CORRECTION2 correct for
;      telluric absorption with the standard DEEP2 routine.  Uses
;      airmass from the FITS header.
;
;      NCORRECT = ncorrect: set this keyword to some threshold (a value
;        of 7, i.e. 7 sigma, seems to work reasonably well; 5 also
;        worked but 10 was too big) to have
;        QE_CORRECTION2 test for and correct high-significance jumps
;        across the chip gap.  Jumps must be at least NCORRECT sigma in
;        significance, as well as being bigger than/smaller than
;        (depending on the sign of the slope)
;        (slope +/- 1.5*NCORRECT*sigma_slope)*delta_lambda, where the slope
;        and robust mean counts are measured using 200 pixels with IVAR
;        greater than 0 at the red end of the blue chip and at the blue
;        end of the red chip.  
;             I.e., we require that there be a significant jump
;        between the blue and red sides, but also that it be much
;        larger than we would expect just due to any slope in the
;        spectrum.  The code errs on the side of not correcting
;        anything.
;             If NCORRECT is not set, no jump correction is applied.
;
;
; PROCEDURES CALLED:
;    MRDFITS
;
; RESTRICTIONS: 
;
;
; EXAMPLE:
;       spec1d_corr=qe_correction2(spec1d, hdr, params, paramsendr, $
;                lambda,throughput_table,NCORRECT=7)
;
;
; MODIFICATION HISTORY:
;
;-APW 9/20/09

; jm14jun02siena - remove this check
;if (max(spec1d.lambda)-min(spec1d.lambda)) gt 3000. $
;   OR mean(spec1d.lambda) lt 7300 OR mean(spec1d.lambda) gt 8300 then begin
;	print,'WARNING: QE_CORRECTION2 will only work for DEEP2-like data'
;        print
;	print,'Returning input array unaltered'
;	return,spec1d
;     endif

if n_elements(header) eq 0 then begin
        print,'You must provide a FITS header'
        print
	print,'Returning input array unaltered'
	return,spec1d
endif

if n_elements(telluric) eq 0 then telluric=0



; remove the a-band atmospheric absorption.
  if telluric then begin
      airmass = sxpar(header, 'AIRMASS')
      remove_telluric, spec1d, airmass;,silent=silent
  endif




;check keyword settings, if empty, reads in the throughput_table,
;params and paramsendr
if n_elements(lambda) eq 0 OR n_elements(throughput_table) eq 0 then begin
;   readcol,  getenv('IDLSPEC1D_DIR')+'/etc/thr_go1200_80_og550.asc',f = 'D,D', $
    readcol,  deep2_path(/auxfiles)+'thr_go1200_80_og550.asc',f = 'D,D', $
        lambda, throughput_table
endif 

if n_elements(params) eq 0 then $
;   params=mrdfits(getenv('IDLSPEC1D_DIR')+'/etc/params.fits')
    params=mrdfits(deep2_path(/auxfiles)+'params.fits')

if n_elements(paramsendr) eq 0 then $
;   paramsendr=mrdfits(getenv('IDLSPEC1D_DIR')+'/etc/paramsendr.fits')
    paramsendr=mrdfits(deep2_path(/auxfiles)+'paramsendr.fits')

num = sxpar(header,'CHIPNO') - 1

npixel=n_elements(spec1d.lambda)

; if array contains full spectrum (not just B or R side):
;   if npixel ge 8192 then begin
     left=indgen(4096)
     right=indgen(4096)+npixel-4096



; divide up wavelengths amongst B and R sides
     xx_b=spec1d.lambda[left]
     xx_r=spec1d.lambda[right]


 
; Ricardo's throughput table seems reasonably reliable when
; tested with stars and coadds from 6300-9300 AA.
whuse=where(lambda gt 6300 and lambda lt 9300)

throughput = interpol(throughput_table[whuse],lambda[whuse],spec1d.lambda)

tmpind = where(spec1d.lambda lt min(lambda[whuse]), ct_add)
if ct_add gt 0 then throughput[tmpind]= throughput[tmpind] > throughput_table[whuse[0]]


; evaluate quadratic fit to correction vs. wavelength given array of
; coefficients, paramsendr
xravg = 8900
yravg = 150

; ramp up correction at the red end of the spectrum, using tabulated
; parameters based on coadded spectra
correctionavg = paramsendr[0] + paramsendr[1]*spec1d.lambda
xavg = (((spec1d.lambda - xravg)/yravg) > 0) < 1
cor2avg = correctionavg*xavg + 1*(1-xavg)
cor2avg = cor2avg > 1

; apply correction for QE differences, based on parameters tabulated
; from comparisons of coadded spectra on different chips
corr_b = params[num,0] + params[num,1]*spec1d.lambda[left] + params[num,2]*spec1d.lambda[left]^2
if num gt 3 then corr_r = right*0+1 else $ ; jm14jun02siena
corr_r = params[num+4,0] + params[num+4,1]*spec1d.lambda[right] + params[num+4,2]*spec1d.lambda[right]^2


; correction was best fit by 1/a quadratic, rather than a quadratic
  corr_b = 1/corr_b
  corr_r = 1/corr_r

  spec1dtmp=spec1d
; apply blue correction
  spec1dtmp.spec[left]=spec1d.spec[left]*corr_b


  spec1dtmp.ivar[left]=spec1d.ivar[left]/(corr_b*corr_b)

; apply red correction
 spec1dtmp.spec[right]=(spec1d.spec[right]*corr_r * cor2avg[right])

 spec1dtmp.ivar[right]=spec1d.ivar[right]/(corr_r*corr_r * cor2avg[right]*cor2avg[right]) 

  spec1d=spec1dtmp

;if necessary, correct for a  jump between the chips
if n_elements(ncorrect) gt 0 then begin
; jm14jun02siena - make sure ivar>0 otherwise linfit() barfs, below
   whl = where(spec1d.ivar[left] gt 0,ct_left)  
  whr = where(spec1d.ivar[right] gt 0,ct_right) 
; whl = where(spec1d.ivar[left] ne 0,ct_left)  
; whr = where(spec1d.ivar[right] ne 0,ct_right) 
  if ncorrect gt 0  and ct_left gt 0 and ct_right gt 0 then begin
   n = ncorrect

;storing the index values where ivar does not equal zero
   whokl = left[whl]
   whokr = right[whr]
   whokall = where(spec1d.ivar ne 0)  

;creating a range of index values around the chip jump
   blueidx = whokl[n_elements(whokl)-321 : n_elements(whokl) - 21]
   redidx = whokr[20 : 319]

;the median counts for the blue and red sides
 
; uses the Hodges-Lehmann mean, a robust mean estimator with better
; performance than the median

   minsig=djsig(spec1d.spec[whokl]) > djsig(spec1d.spec[whokr])

   ctsb = hlmean(spec1d.spec[blueidx]) 
   sigmab = stdev(spec1d.spec[blueidx]) > minsig
   ctsr = hlmean(spec1d.spec[redidx])
   sigmar = stdev(spec1d.spec[redidx])  > minsig

;the difference in counts from red to blue
   delta = ctsr - ctsb

   sigma_medianb = sigmab / (sqrt((.9) * 300))
   sigma_medianr = sigmar / (sqrt((.9) * 300))

; significance of the difference
   delta_sig = abs(delta) / sqrt((sigma_medianb)^2 + (sigma_medianr)^2)

;first flag : is the difference significant?
   if delta_sig gt n then tag1 = 1 else tag1 = 0


;the change in the average wavlengths of the two regions
   lambda_blue= mean(spec1d.lambda[blueidx])
   lambda_red= mean(spec1d.lambda[redidx])
   delta_lambda = lambda_red - lambda_blue

;find the slope of each region
; jm14jun02siena - use linfit()    
   paramsb = linfit(spec1d.lambda[blueidx], spec1d.spec[blueidx],sigma=sigmab,$
          measure=1/(sqrt(spec1d.ivar[blueidx])))
   paramsr = linfit(spec1d.lambda[redidx],spec1d.spec[redidx],sigma=sigmar, $
          measure=1/(sqrt(spec1d.ivar[redidx])))
;   paramsb = svdfit(spec1d.lambda[blueidx], spec1d.spec[blueidx],2,sigma=sigmab,$
;          measure=1/(sqrt(spec1d.ivar[blueidx])),yfit=yfitb)
;   paramsr = svdfit(spec1d.lambda[redidx],spec1d.spec[redidx] , 2,sigma=sigmar, $
;          measure=1/(sqrt(spec1d.ivar[redidx])))
   
; figure out the sign of the slope
   signb=paramsb[1]/abs(paramsb[1])
   signr=paramsr[1]/abs(paramsr[1])

;flag 2, if the change in average counts from blue to red is greater
;than the predicted change in counts from the slope, augmented by
;1.5*NCORRECT sigma 


if delta gt 0 then begin
    flagb = delta gt (paramsb[1]+1.5*N*sigmab[1])*delta_lambda
    flagr = delta gt (paramsr[1]+1.5*N*sigmar[1])*delta_lambda 
endif else begin
    flagb = delta lt (paramsb[1]-1.5*N*sigmab[1])*delta_lambda
    flagr = delta lt (paramsr[1]-1.5*N*sigmar[1])*delta_lambda 
endelse 

 ;  if signb gt 0 then flagb = delta gt (paramsb[1]+1.*N*sigmab[1])*delta_lambda $
 ;       else flagb = delta lt (paramsb[1]-1.*N*sigmab[1])*delta_lambda

 ;  if signr gt 0 then flagr = delta gt (paramsr[1]+1.*N*sigmar[1])*delta_lambda $
 ;       else flagr = delta lt (paramsr[1]-1.*N*sigmar[1]) * delta_lambda

; only correct if the change is significant using the slope from
; either the blue or red side
   tag2 = flagb AND flagr

; if the change is significant, correct the blue to the red side
   if (tag1 and tag2) then spec1d.spec[left]=spec1d.spec[left]*ctsr/ctsb
   if (tag1 and tag2) then spec1d.ivar[left]=spec1d.ivar[left]/(ctsr/ctsb)^2
   if (tag1 AND tag2) then print,'correcting chip jump of ',$
    string(delta_sig,form='(F6.2)'),' sigma significance,multiplying by ',string(ctsr/ctsb,format='(F6.3)')

  endif
endif
  spec1dtmp=spec1d


  spec1dtmp.spec[left]=spec1d.spec[left]/throughput[left]
  spec1dtmp.ivar[left]=spec1d.ivar[left]*throughput[left]^2

; apply red correction
 spec1dtmp.spec[right]=spec1d.spec[right] / throughput[right]

 spec1dtmp.ivar[right]=spec1d.ivar[right] * throughput[right]^2

return,spec1dtmp

end
