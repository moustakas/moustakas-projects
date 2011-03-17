pro ages_inspect_specfit, pass
; jm09mar06nyu

    npass = n_elements(pass)
    if (npass eq 0L) then pass = 101
    pass = string(pass,format='(I0)')

    aa1 = read_ages(/ancillary)
    for ipass = 0l, npass-1l do begin
    
       ww = where(pass[ipass] eq aa1.pass)
       aa = aa1[ww]

       ages_display_spectrum, aa, /plotobswave, plottype=2
       
;      ss = read_ages_specfit(aa.galaxy)

    endfor

stop
    
return
end
    
