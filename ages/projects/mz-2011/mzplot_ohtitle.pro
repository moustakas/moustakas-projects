function mzplot_ohtitle, t04=t04, m91=m91, kk04=kk04, fluxcor=fluxcor, mean=mean
    calib = ''
    if keyword_set(fluxcor) then calib = 'cor'
    if keyword_set(t04) then calib = 'T04'
    if keyword_set(m91) then calib = 'M91'
    if keyword_set(kk04) then calib = 'KK04'

    if (keyword_set(t04) and keyword_set(fluxcor)) then calib = 'T04,cor'
    if (keyword_set(m91) and keyword_set(fluxcor)) then calib = 'M91,cor'
    if (keyword_set(kk04) and keyword_set(fluxcor)) then calib = 'KK04,cor'

    if (keyword_set(t04) and keyword_set(fluxcor) eq 0) then calib = 'T04,EW'
    if (keyword_set(m91) and keyword_set(fluxcor) eq 0) then calib = 'M91,EW'
    if (keyword_set(kk04) and keyword_set(fluxcor) eq 0) then calib = 'KK04,EW'

    title = textoidl('12 + log (O/H)_{'+calib+'}')
    if keyword_set(mean) then title = '<'+title+'>'
   
return, title
end
    
