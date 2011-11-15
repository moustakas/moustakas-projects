function mzplot_sfrmtitle, avg=avg
    H0 = string(100*mz_h100(),format='(I0)')
    if keyword_set(avg) then return, textoidl('<log (\psi / M)> (Gyr^{-1})')
    return, textoidl('log (\psi / M) (Gyr^{-1})')
end
