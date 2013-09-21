function read_assessment, file, questions=qq, key=key
; jm13sep04siena - read the CSV-formated GoogleDocs spreadsheet;
; assumes just 30 questions on the FCI!

    data1 = djs_readlines(file,nhead=1,head=header)
    nstudent = n_elements(data1)
    nkey = n_elements(key)
    
; deal with blank data and also commas in the long-answers; this is
; stupidly hard-coded to probably three blank answers
    data2 = repstr(repstr(repstr(data1,', ',' '),',,',',...,'),',,',',...,')
    ncol = n_elements(strsplit(data2[0],',',/extract))

; parse the questions
    qq = strsplit(header,',',/extract)

; and now the data
    data = {timestamp: '', email: ''}
    for ii = 2, ncol-1 do data = create_struct(data,'Q'+$
      string(ii-2+1,format='(I2.2)'),'')
    data = replicate(data,nstudent)

    for jj = 0, nstudent-1 do begin
       if strmid(data2[jj],0,strlen(data2[jj])-1,/reverse) eq ',' then $ ; last column is empty
         data2[jj] = data2[jj]+'...'
       str = strsplit(data2[jj],',',/extract)
       if nkey ne 0 then begin
          rightwrong1 = str[2:31] eq key
          if jj eq 0 then rightwrong = rightwrong1 else $
            rightwrong = [[rightwrong],[rightwrong1]]
       endif
       for ii = 0, ncol-1 do data[jj].(ii) = str[ii]
    endfor

; if given the key, add the final score and the distribution of
; right/wrong answers 
    if nkey ne 0 then begin
       data = struct_addtags(data,replicate({score: 0.0, rightwrong: bytarr(nkey)},nstudent))
       data.rightwrong = rightwrong
       data.score = total(rightwrong,1)/float(nkey)
    endif

return, data
end

