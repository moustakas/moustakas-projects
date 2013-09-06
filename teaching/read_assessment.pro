function read_assessment, file, questions=qq
; jm13sep04siena - read the CSV-formated GoogleDocs spreadsheet 

    data1 = djs_readlines(file,nhead=1,head=header)
    nstudent = n_elements(data1)

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
       for ii = 0, ncol-1 do data[jj].(ii) = str[ii]
    endfor
       
return, data
end

