pro ks_function, mags

        common trgb_ks, alph, bet, trgb, width

        x = mags-trgb
        model = 10.^(alph*x)*(x ge 0.) + 10.^(-width+bet*x)*(x lt 0)

return, model        
end
