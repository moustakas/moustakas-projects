function sfrm_masstitle, mstar=mstar
    if keyword_set(mstar) then return, $
      textoidl('log (!8M^{*}!6h_{70}^2/!8M'+sunsymbol()+'!6)')
    return, textoidl('log_{10}(!8M!6h_{70}^2/!8M'+sunsymbol()+'!6)')
;   return, textoidl('log_{10}(!8M!6/!8M'+$
;     sunsymbol()+'!6) + log_{10}(h_{70}^{2})')
end
    
