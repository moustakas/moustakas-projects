pro astr010_s18_grades, alldata, test=test, sendit=sendit, final=final, drop=drop
; jm18jan29siena - parse the grades for this class 

    path = getenv('TEACHING_DIR')+'/010-S18/grades/'
    
    date = '18mar09' ; update this
    semester = 'Spring 2018'
    class = 'ASTRO010 - Introductory Astronomy'

; read the grade spreadsheet downloaded from GoogleDocs
    gradefile = path+'astr010_grades_'+date+'.csv'
    data = read_gradefile(gradefile,unique_assignments=assign)

; specify the complete list of *possible* assignments and their
; relative weights
;    allassign = ['Homework','Module Papers','Module Exams','Final Exam']
;    weight = [0.20,0.30,0.30,0.20]
    allassign = ['Homework','Module I Paper','Module II Paper','Module III Paper',$
      'Module I Exam','Module II Exam','Module III Exam','Final Exam']
    weight = [0.20,0.10,0.10,0.10,0.10,0.10,0.10,0.20]
    if keyword_set(drop) then begin
       droplowest = [1,0,0,0,0,0,0,0]
    endif else begin
       droplowest = [0,0,0,0,0,0,0,0]
    endelse

    process_grades, data, assign=assign, allassign=allassign, $
      weight=weight, class=class, semester=semester, test=test, $
      sendit=sendit, alldata=alldata, final=final, droplowest=droplowest
    srt = sort(alldata.current_grade)
    struct_print, alldata[srt]

return
end
