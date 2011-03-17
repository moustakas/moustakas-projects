pro ir_plots, ir=ir
; jm05aug15uofa

    if (n_elements(ir) eq 0L) then ir = read_irgalaxies(/templates_04)
    nir = n_elements(ir)

    for i = 30, nir-1L do begin
;   for i = 0L, nir-1L do begin

       yesmass = where(ir[i].continuum_mass_fraction gt 0.0)

       age = (ir[i].template_age[yesmass]>1E5)/1E6 ; NOTE!
       fraction = 100*(ir[i].continuum_mass_fraction[yesmass]>1E-3)

       title = strtrim(ir[i].nice_galaxy,2)+' [log (M/M'+sunsymbol()+')='+strtrim(string(alog10(ir[i].continuum_total_mass),$
         format='(F12.2)'),2)+']'
       
       plot, age, fraction, ps=-4, /ylog, thick=2, charsize=1.6, $
         charthick=2, xsty=3, ysty=3, title=title, $
         yrange=100*[1E-3,1.0], xrange=[1E5,1.2E10]/1E6, /xlog, $
         ytitle='Mass Fraction [%]', xtitle='Time [Myr Ago]'
       cc = get_kbrd(1)

    endfor
    
    
return
end
    
