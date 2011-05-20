function mzplot_sfrmtitle
    H0 = string(100*mz_h100(),format='(I0)')
    return, textoidl('log (\psi / M) (Gyr^{-1})')
end
