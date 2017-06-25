pro astr392_s17_grades, alldata, test=test, sendit=sendit, final=final
; jm15jan31siena - parse the grades for this class 

    path = getenv('TEACHING_DIR')+'/392-S17/grades/'
    
    date = '17may08' ; update this
    semester = 'Spring 2017'
    class = 'Astronomy 392 - Principles of Astrophysics II'

; read the grade spreadsheet downloaded from GoogleDocs
    gradefile = path+'astr392_grades_'+date+'.csv'
    data = read_gradefile(gradefile,unique_assignments=assign)

; specify the complete list of *possible* assignments and their
; relative weights
    allassign = ['Midterm 1 Exam','Midterm 2 Exam',$
      'Midterm 2 Computational','Astrophysics Talk','Final Problem Set']
    weight = [0.25,0.2125,0.0375,0.25,0.25]
    droplowest = [0,0,0,0,0]

; include this factor in the final grades --    
;   data.final_problem_set *= 1.07
    
;   keep = where(strmatch(data.last_name,'*Young*'))
;   keep = where(data.final_exam gt 0.0)
;   data = data[keep]

    process_grades, data, assign=assign, allassign=allassign, $
      weight=weight, class=class, semester=semester, test=test, $
      sendit=sendit, alldata=alldata, final=final, droplowest=droplowest
    srt = sort(alldata.current_grade)
    struct_print, alldata[srt]

return
end
