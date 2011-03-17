function sfrm_phititle, phistar=phistar
    if keyword_set(phistar) then return, $
      textoidl('\Phi^{*}(!8M!6) (10^{-3} h_{70}^3 Mpc^{-3} dex^{-1})') 
    return, textoidl('\Phi(!8M!6) (h_{70}^3 Mpc^{-3} dex^{-1})') 
end
    

