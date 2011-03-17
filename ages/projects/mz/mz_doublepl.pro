function mz_doublepl, logmass, params
; return a smoothly connected double power-law
    logohstar = params[0]
    logmstar = params[1]
    alpha1 = params[2]
    alpha2 = params[3]

    ohstar = 10D^(logohstar-12.0)
    model = ohstar / ((10D^(logmass-logmstar))^alpha1 + (10D^(logmass-logmstar))^alpha2)
return, 12+alog10(model)
end

