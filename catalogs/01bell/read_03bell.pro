function read_03bell, johnson=johnson
; jm07feb06nyu - excised from READ_01BELL

; read Table 7 of Bell et al. 2003

    root = '01bell'
    
    path = getenv('CATALOGS_DIR')+'/'+root+'/'
    if keyword_set(johnson) then $
      bell03 = rsex(path+'03bell_table7_johnson.dat') else $
        bell03 = rsex(path+'03bell_table7.dat')
    
return, bell03
end
