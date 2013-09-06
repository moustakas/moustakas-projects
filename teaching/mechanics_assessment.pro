;function fci_key
;; correct answers
;    key = [$
;1 B 
;2 D 
;3 E 
;4 C 
;5 A 
;6 C 
;7 C 
;8 D 
;9 A 
;10 E 
;11 E 
;12 C 
;13 B 
;14 B 
;15 E 
;16 A 
;17 D 
;18 B 
;19 C 
;20 C 
;21 A 
;22 B 
;23 D 
;24 A 
;25 A 
;26 E

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
    path = getenv('TEACHING_DIR')+'/assessment/'+year+'-mechanics/'
    testfile = file_search(path+'*.csv',count=ntest)

    bin = 0.5
    
    for ii = 0, ntest-1 do begin
       data = read_assessment(testfile[ii],questions=qq) ; read it
       nqq = n_elements(qq)

       for jj = 2, nqq-1 do begin
          im_plothist, remap_abcde(data.(jj))-bin, bin=bin, edge=1.0, $
            xrange=[0.5,5.5], xsty=1, /fill, xtickname=['A','B','C','D','E'], $
            ytitle='Number of Students'
          
          cc = get_kbrd(1)
       endfor
          
    endfor
       

stop
    
return
end
