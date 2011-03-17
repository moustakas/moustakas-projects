; jm00nov15uofa
; take a TRGB magnitude, determine the distance modulus, and output
; the distance and the error in the distance
pro dmod2d, tip, tiperr, a_i, dmod, sigdmod, d, sigd

	itip = -4.06     ; absolute I_TRGB 
        sigitip = 0.07   ; random error

        dmod = tip-itip+a_i
        sigdmod = sqrt(tiperr^2+sigitip^2)

        d = 10.0^((dmod-25.0)/5.0)     ; Mpc
        sigd = alog(10)*d*sigdmod/5.0  ; error

return
end
