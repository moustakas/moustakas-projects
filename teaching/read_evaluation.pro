function read_evaluation, file, questions=qq
; jm12oct19siena - read a mid-semester evaluation as a CSV-formated
;   GoogleDocs spreadsheet

    data1 = djs_readlines(file,nhead=1,head=header)
    nstudent = n_elements(data1)

; deal with blank data and also commas in the long-answers; this is
; stupidly hard-coded to probably three blank answers
    data2 = repstr(repstr(repstr(data1,', ',' '),',,',',...,'),',,',',...,')
    ncol = n_elements(strsplit(data2[0],',',/extract))

; parse the questions
    qq = strsplit(header,',',/extract)

; and now the data
    data = strarr(ncol,nstudent)
    for ii = 0, nstudent-1 do begin
       if strmid(data2[ii],0,strlen(data2[ii])-1,/reverse) eq ',' then $ ; last column is empty
         data2[ii] = data2[ii]+'...'
       str = strsplit(data2[ii],',',/extract)
;      splog, ii & help, qq, str & print
       data[*,ii] = str
    endfor

return, data
end

