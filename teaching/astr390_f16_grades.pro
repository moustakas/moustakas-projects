pro astr390_f16_grades, alldata, test=test, sendit=sendit, final=final, drop=drop
; jm15jan31siena - parse the grades for this class 

    path = getenv('TEACHING_DIR')+'/390-F16/grades/'
    
    date = '16dec20' ; update this
    semester = 'Fall 2016'
    class = 'ASTR390 - Principles of Astrophysics I'

; read the grade spreadsheet downloaded from GoogleDocs
    gradefile = path+'astr390_grades_'+date+'.csv'
    data = read_gradefile(gradefile,unique_assignments=assign)

; specify the complete list of *possible* assignments and their
; relative weights
    allassign = ['Homework','Quizzes','Talk','Final Problem Set']
    weight = [0.30,0.40,0.15,0.15]
    if keyword_set(drop) then begin
       droplowest = [1,1,0,0]
    endif else begin
       droplowest = [0,0,0,0]
    endelse

;   keep = where(strmatch(data.first_name,'*Fred*'))
;   data = data[keep]

;    if keyword_set(final) eq 0 then begin
;       data.final_problem_set *= 1.1
;       data.quizzes *= 1.05
;    endif

    process_grades, data, assign=assign, allassign=allassign, $
      weight=weight, class=class, semester=semester, test=test, $
      sendit=sendit, alldata=alldata, final=final, droplowest=droplowest
    srt = sort(alldata.current_grade)
    struct_print, alldata[srt]

return
end
