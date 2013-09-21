function get_title, qq
    words = strsplit(qq,' ',/extract)
    nw = n_elements(words)
    if nw lt 5 then tit = qq else $
      tit = strjoin([words[0:nw/2-1],'!c',words[nw/2:nw-1]],' ')
return, tit    
end

function remap_abcde, str
    out = str
    out = repstr(out,'A','1')
    out = repstr(out,'B','2')
    out = repstr(out,'C','3')
    out = repstr(out,'D','4')
    out = repstr(out,'E','5')
return, float(out)
end

pro mechanics_assessment
; jm13sep04siena

    year = ['2013']
    rootpath = getenv('TEACHING_DIR')+'/assessment/mechanics/'
    filepath = rootpath+year+'/'
    
; read the answer key and then all the data
    readcol, rootpath+'mechanics-key.txt', key, format='A', /silent
    nkey = n_elements(key) ; =30 questions
    
    testfile = file_search(filepath+'*.csv',count=ntest)

; read all the data
    for ii = 0, ntest-1 do begin
       data1 = read_assessment(testfile[ii],key=key,questions=qq)
       if ii eq 0 then data = data1 else data = [data,data1]
    endfor
    nall = n_elements(data)
    
    

stop

    bin = 0.5
    for ii = 0, ntest-1 do begin
       data = read_assessment(testfile[ii],key=key,$ ; read it
         questions=qq,rightwrong=rightwrong) 
       
       nqq = n_elements(qq)

       for jj = 2, 2+nkey-1 do begin
          im_plothist, remap_abcde(data.(jj))-bin, bin=bin, edge=1.0, $
            xrange=[0.5,5.5], xsty=1, /fill, xtickname=['A','B','C','D','E'], $
            ytitle='Number of Students'
          
          cc = get_kbrd(1)
       endfor
          
    endfor


return
end
