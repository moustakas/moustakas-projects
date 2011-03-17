pro compare_stasinska_garnett
; jm04jul13uofa
; what is the systematic difference between T(O+) temperatures
; predicted by Stasinska (1990) versus Garnett (1992)?

    tmax = 22000.0
    tmin = 4000.0
    dt = 100.0
    
    toiii = findgen((tmax-tmin)/dt+1)*dt+tmin

    toii_g92 = 0.7*toiii + 3000.0
    toii_s90 = 20000.0/(10000.0/toiii + 0.8)

    presid = 100.0*(toii_g92-toii_s90)/toii_g92
    
    plot, toiii, presid, line=2, xrange=[tmin,tmax], xsty=3, ysty=3, $
      xthick=2.0, ythick=2.0, thick=2.0, charthick=2.0, charsize=2.0, $
      xtitle='T([O III]) [K]', ytitle='Percent Residuals [(G92 - S90) / G92]', $
      xtickname=['5000','10000','15000','20000']
    oplot, !x.crange, [0,0], line=0, thick=2.0
    
stop    

return
end
    
