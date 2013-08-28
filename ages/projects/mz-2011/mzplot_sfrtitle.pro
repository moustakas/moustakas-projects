function mzplot_sfrtitle
    H0 = string(100*mz_h100(),format='(I0)')
    return, textoidl('log (\psi / M_{'+sunsymbol()+'} yr^{-1})')
end
