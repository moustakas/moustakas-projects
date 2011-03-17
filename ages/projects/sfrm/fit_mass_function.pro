function fit_mass_function, mass, weight, binsize=binsize, $
  binmass=binmass, phi=phi, errphi=phierr, parinfo=parinfo, $
  nofit=nofit, plotlimit=plotlimit, masslimit=masslimit
; fit the stellar mass function

; round the bin centers to two decimal points
    histmin = fix((min(mass)+binsize/2.0)*100.0)/100.0-binsize/2.0

; build the MF for the plot, setting the weight to zero for objects
; below PLOTLIMIT
    phi = im_hist1d(mass,weight*(mass gt plotlimit),$
      binsize=binsize,obin=binmass,binedge=0,$
      h_err=phierr,histmin=histmin)
    
; fit the mass function; remake the mass function, setting the weight
; to zero for objects less massive than MASSLIMIT
    if keyword_set(nofit) then fit = 0.0 else begin
       fitphi1 = im_hist1d(mass,weight*(mass gt masslimit),$
         binsize=binsize,obin=fitmass1,binedge=0,$
         h_err=fitphierr1,histmin=histmin)
       good = where(fitphi1 ne 0.0)
       fitmass = fitmass1[good]
       fitphi = fitphi1[good]
       fitphierr = fitphierr1[good]

       mf_fit_schechter_plus, 10^fitmass, fitphi, $
         fitphierr, fit, parinfo=parinfo

;      ploterror, binmass, phi, phierr, ps=10, xsty=3, ysty=3, $
;        yrange=[1E-5,1E-1], /ylog
;      oploterror, fitmass, fitphi, fitphierr, ps=10, $
;        color=djs_icolor('red'), errcolor=djs_icolor('red')
;      djs_oplot, fitmass, mf_schechter(10^fitmass,fit), color='yellow'
;      djs_oplot, fitmass, mf_schechter(10^fitmass,parinfo[0].value,$
;        parinfo[1].value,parinfo[2].value), color='green'
    endelse
return, fit
end

