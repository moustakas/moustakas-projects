function read_01bell
; jm05may17uofa

; read Table 3 of Bell & de Jong 2001

    root = '01bell'
    
    path = getenv('CATALOGS_DIR')+'/'+root+'/'
    bell01 = im_read_fmr(path+root+'_table3.dat')

    keep = lindgen(6)*7+5
    bell01 = bell01[keep]
;   niceprint, bell01.color, bell01.model
    
return, bell01
end
