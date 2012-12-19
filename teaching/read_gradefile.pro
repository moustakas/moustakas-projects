function read_gradefile, file, unique_assignments=uassign
; jm12oct05siena - read a CSV-formated GoogleDocs spreadsheet
;   containing the scores for each assignment

; read and parse the grade matrix and the four-line header
    nhead = 4
    data1 = djs_readlines(file,nhead=nhead,head=header1)
    data1 = repstr(data1,'--','0')
    nstudent = n_elements(data1)

    ncol = n_elements(strsplit(data1[0],',',/extract))
    data = strarr(ncol,nstudent)
    header = strarr(ncol,nhead)

    for ii = 0, nstudent-1 do data[*,ii] = strsplit(data1[ii],',',/extract)
    for ii = 0, nhead-1 do header[*,ii] = strsplit(header1[ii],',',/extract)

; assume that the first FIVE columns contain the student data (INFO),
; and the subsequent columns contain the assignment types (ASSIGN)
    ninfo = 6
    nassign = ncol-ninfo
    sdata = data[0:ninfo-1,*]
    data = data[ninfo:ncol-1,*]
    sheader = header[0:ninfo-1,0] ; the rest of the lines are blank
    header = header[ninfo:ncol-1,*]

; get the full and unique list of assignments
    assign = strtrim(header[*,0],2)
    uassign = assign[uniq(assign,sort(assign))]
    nuassign = n_elements(uassign)

; build the output structure    
    out = create_struct(sheader[0],' ') ; should be "last name"
    for ii = 1, ninfo-1 do out = create_struct(out,sheader[ii],'')
    for ii = 0, nuassign-1 do begin
       match = where(uassign[ii] eq assign,nmatch)
       nmatch = strtrim(nmatch,2)

       tags = [uassign[ii],uassign[ii]+'_details',uassign[ii]+'_date',uassign[ii]+'_points']
       if nmatch eq 1 then $
         types = ['F','A','A','F'] else $
           types = ['F('+nmatch+')','A('+nmatch+')','A('+nmatch+')','F('+nmatch+')']

       create_struct, junk, '', tags, types
       out = create_struct(out,junk)

       indx = tag_indx(out,repstr(uassign[ii]+'_details',' ','_')) ; second row
       out.(indx) = header[match,1]

       indx = tag_indx(out,repstr(uassign[ii]+'_date',' ','_')) ; third row
       out.(indx) = header[match,2]

       indx = tag_indx(out,repstr(uassign[ii]+'_points',' ','_')) ; fourth row
       out.(indx) = header[match,3]
    endfor
    out = replicate(out,nstudent)

; parse the student info
    for ii = 0, ninfo-1 do out.(ii) = reform(sdata[ii,*])

; now parse the points earned for each assignment
    for ii = 0, nuassign-1 do begin
       match = where(uassign[ii] eq assign,nmatch)
       indx = tag_indx(out,repstr(uassign[ii],' ','_'))
       out.(indx) = reform(data[match,*])
    endfor

return, out
end

