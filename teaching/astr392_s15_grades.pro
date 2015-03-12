pro astr392_s15_grades, alldata, test=test, sendit=sendit, final=final
; jm15jan31siena - parse the grades for this class 

    path = getenv('TEACHING_DIR')+'/392-S15/grades/'
    
    date = '15mar08' ; update this
    semester = 'Spring 2015'
    class = 'Astronomy 392 - Principles of Astrophysics II'

; read the grade spreadsheet downloaded from GoogleDocs
    gradefile = path+'astr392_grades_'+date+'.csv'
    data = read_gradefile(gradefile,unique_assignments=assign)

; specify the complete list of *possible* assignments and their
; relative weights
    allassign = ['Homework','Quizzes','Final Problem Set']
    weight = [0.40,0.40,0.20]
    droplowest = [0,1,0]

;   keep = where(strmatch(data.last_name,'*Mej*'))
;   keep = where(data.final_exam gt 0.0)
;   data = data[keep]

    process_grades, data, assign=assign, allassign=allassign, $
      weight=weight, class=class, semester=semester, test=test, $
      sendit=sendit, alldata=alldata, final=final, droplowest=droplowest
    srt = sort(alldata.current_grade)
    struct_print, alldata[srt]

return
end
