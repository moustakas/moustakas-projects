function mz_closedbox, logmass, params
; see the discussion section in the paper

    mratio = params[1]/10D^logmass
;   mratio = params[1]/(10^(logmass-9D))
; good    
;   model = params[0] - mratio^(-params[2])
    model = params[0] - alog10(1+mratio^params[2])

; tests
;   model = 12.0+alog10(params[0]/(1+mratio^params[2]))
;   model = params[0] - mratio^(-params[2])
;   model = params[0] - params[2]*alog10(mratio) + params[3]*mratio^(-params[2])
;   model = params[0] + alog10(mratio^(-params[2])) + mratio^(-params[2])
;   model = params[0]*mratio^(-params[2])*exp(-mratio)^(-params[3])
;   model = params[0] + params[1]*(10D^(params[2]-logmass))^params[3] + $
;     params[3]*alog10(10D^(params[2]-logmass))
;   model = params[0] - params[3]*((10D^(logmass-params[1]))^(-params[2]))
return, model
end

