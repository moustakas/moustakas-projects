function cosmicimf_rhotitle, sfr=sfr
    if keyword_set(sfr) then return, textoidl('log (\rho_{SFR}/h_{70} M_{\odot} yr^{-1} Mpc^{-3})')
    return, textoidl('log (\rho_{*}/h_{70} M_{\odot} Mpc^{-3})')
end
