function rykoff_mass, lambda
; relation between richness and M200 taken from Appendix B in
; Rykoff+12 
;   return, alog10(1D14*exp(1.72+1.08*alog(lambda/60.0))) ; B4 (M200m)
    return, alog10(1D14*exp(1.14+1.04*alog(lambda/60.0))) ; B6 (M500c)
end

