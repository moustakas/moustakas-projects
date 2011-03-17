function cosmicimf_phititle, phistar=phistar, log=log
    if keyword_set(log) then begin
       if keyword_set(phistar) then $
         title = textoidl('log (\Phi^{*}/h_{70}^3 Mpc^{-3} dex^{-1})'); else $
;          title = textoidl('\Phi(!8M!6) (h_{70}^3 Mpc^{-3} dex^{-1})')
    endif else begin
       if keyword_set(phistar) then $
         title = textoidl('\Phi^{*}(M) (10^{-3} h_{70}^3 Mpc^{-3} dex^{-1})') else $
;        title = textoidl('\Phi^{*}(!8M!6) (10^{-3} h_{70}^3 Mpc^{-3} dex^{-1})') else $
           title = textoidl('\Phi(M) (h_{70}^3 Mpc^{-3} dex^{-1})')
;          title = textoidl('\Phi(!8M!6) (h_{70}^3 Mpc^{-3} dex^{-1})')
    endelse

return, title
end
    

