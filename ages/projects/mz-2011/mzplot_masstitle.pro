function mzplot_masstitle, mstar=mstar
    H0 = string(100*mz_h100(),format='(I0)')
    if keyword_set(mstar) then return, $
      textoidl('log (M^{*} / h_{'+H0+'}^{-2} M'+sunsymbol()+')')
    return, textoidl('log (M / h_{'+H0+'}^{-2} M'+sunsymbol()+')')

;   if keyword_set(mstar) then return, $
;     textoidl('log (!8M^{*}!6h_{70}^2/!8M'+sunsymbol()+'!6)')
;   return, textoidl('log (!8M!6h_{70}^2/!8M'+sunsymbol()+'!6)')
;   return, textoidl('log_{10}(!8M!6/!8M'+$
;     sunsymbol()+'!6) + log_{10}(h_{70}^{2})')
end
