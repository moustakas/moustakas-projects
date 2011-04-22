pro clash_merge_cat
; jm11apr11ucsd - merge the catalogs

    allcat = file_search('*.cat',count=ncat)
    for ii = 0, ncat-1 do begin
       cat1 = rsex(allcat[ii])
       filt = strmid(file_basename(allcat[ii]),0,strpos(file_basename(allcat[ii]),'_'))
       
stop
    endfor
    
    
return
end
    
