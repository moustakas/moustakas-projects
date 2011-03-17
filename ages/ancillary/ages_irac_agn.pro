function ages_irac_agn, ages, debug=debug
; jm10mar16ucsd - select AGN using the Stern diagram

    v2ab = k_vega2ab(filterlist=irac_filterlist(),/kurucz,/silent)*0.0 ; already Vega

    isagn = intarr(n_elements(ages))-1 ; 
    good = where($
      (ages.ch1_mag_aper_04 gt 0.0) and (ages.ch1_mag_aper_04 lt 90.0) and $
      (ages.ch2_mag_aper_04 gt 0.0) and (ages.ch2_mag_aper_04 lt 90.0) and $
      (ages.ch3_mag_aper_04 gt 0.0) and (ages.ch3_mag_aper_04 lt 90.0) and $
      (ages.ch4_mag_aper_04 gt 0.0) and (ages.ch4_mag_aper_04 lt 90.0),ngood)

    ch12 = (ages[good].ch1_mag_aper_04-ages[good].ch2_mag_aper_04)-(v2ab[0]-v2ab[1])
    ch34 = (ages[good].ch3_mag_aper_04-ages[good].ch4_mag_aper_04)-(v2ab[2]-v2ab[3])

    isagn[good] = (ch34 gt 0.6) and $
      (ch12 gt (0.2*ch34+0.18)) and $
      (ch12 gt (2.5*ch34-3.5))

    if keyword_set(debug) then begin
       djs_plot, ch34, ch12, psym=6, sym=0.2, xsty=3, ysty=3, $
         xrange=[-0.5,3.5], yrange=[-0.3,1.5]
       ch34axis = im_array(-1.0,5.0,0.02)
       djs_oplot, ch34axis, poly(ch34axis,[0.18,0.2])
       djs_oplot, ch34axis, poly(ch34axis,[-3.5,2.5])
       djs_oplot, 0.6*[1,1], !y.crange
       agn = where(isagn[good],nagn)
       if (nagn ne 0L) then djs_oplot, ch34[agn], $
         ch12[agn], psym=6, sym=0.2, color='red'
    endif

return, isagn
end

