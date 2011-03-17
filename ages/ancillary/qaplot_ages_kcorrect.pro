pro qaplot_ages_kcorrect, info, psfile=psfile, vname=vname, clobber=clobber
; jm10feb10ucsd - wrapper on KCORRECT_QAPLOT; INFO should be something
;   analogous to the structure written out by AGES_KCORRECT

    if (n_elements(vname) eq 0) then vname = 'default.nolines'
    ngal = n_elements(info)
    if (ngal eq 0L) then begin
       doc_library, 'kcorrect_qaplot'
       return
    endif

    select = ['filterlist','bestmaggies','maggies','ivarmaggies','coeffs','z','chi2','mass']
    newtags = ['in_filterlist','bestmaggies','maggies','ivarmaggies','coeffs','z','chi2','mass']
    info1 = im_struct_trimtags(info,select=select,newtags=newtags)
    
    kcorrect_qaplot, info1, psfile=psfile, vname=vname, clobber=clobber

return
end
    
