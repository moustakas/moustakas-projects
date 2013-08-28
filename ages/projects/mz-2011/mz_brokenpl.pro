function mz_brokenpl, logmass, params
; return a broken power-law
    logohstar = params[0]
    logmstar = params[1]
    alpha1 = params[2]
    alpha2 = params[3]

    ohstar = 10D^(logohstar-12.0)
    model = logmass*0.0
    below = where(logmass lt logmstar,nbelow,comp=above,ncomp=nabove)
    if (nbelow gt 0L) then model[below] = ohstar*(10D^(logmass[below]-logmstar))^alpha1
    if (nabove gt 0L) then model[above] = ohstar*(10D^(logmass[above]-logmstar))^alpha2
return, 12+alog10(model)
end

