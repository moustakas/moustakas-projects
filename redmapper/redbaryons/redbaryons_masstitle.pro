function redbaryons_masstitle, mstar=mstar, simple=simple, $
  cumuphi=cumuphi
    
    H0 = string(100*mf_h100(),format='(I0)')

;   title = 'log (M / M'+sunsymbol()+')'
    title = 'log (M_{*} / h_{'+H0+'}^{-2} M'+sunsymbol()+')'
    if keyword_set(mstar) then title = $
      'log (M^{*} / h_{'+H0+'}^{-2} M'+sunsymbol()+')'
    if keyword_set(simple) then title = 'log (M/M'+sunsymbol()+')'
    if keyword_set(cumuphi) then title = 'log (M_{c} / M'+sunsymbol()+')'
;   if keyword_set(cumuphi) then title = 'log (M_{c} / h_{'+H0+'}^{-2} M'+sunsymbol()+')'
    
;   if keyword_set(mstar) then return, $
;     textoidl('log (M^{*}h_{'+H0+'}^2/M'+sunsymbol()+')')
;   return, textoidl('log (Mh_{'+H0+'}^2/M'+sunsymbol()+')')

;   if keyword_set(mstar) then return, $
;     textoidl('log (!8M^{*}!6h_{'+H0+'}^2/!8M'+sunsymbol()+'!6)')
;   return, textoidl('log (!8M!6h_{'+H0+'}^2/!8M'+sunsymbol()+'!6)')
;   return, textoidl('log_{10}(!8M!6/!8M'+$
;     sunsymbol()+'!6) + log_{10}(h_{70}^{2})')

return, textoidl(title)
end
    
