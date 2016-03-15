pro phys140_s16_grades, alldata, test=test, sendit=sendit, final=final
; jm15jan27siena - parse the grades for this class 

    path = getenv('TEACHING_DIR')+'/140-S16/grades/'
    
    date = '16mar14' ; update this
    semester = 'Spring 2016'
    class = 'Physics 140 - General Physics II'

; read the grade spreadsheet downloaded from GoogleDocs
    gradefile = path+'phys140_grades_'+date+'.csv'
    data = read_gradefile(gradefile,unique_assignments=assign)

; specify the complete list of *possible* assignments and their
; relative weights
    allassign = ['Homework','Quizzes','Lab','Final Exam']
    weight = [0.25,0.35,0.15,0.25]
    droplowest = [1,1,0,0]

;   keep = where(strmatch(data.last_name,'*Nap*'))
;   keep = where(data.final_exam gt 0.0)
;   data = data[keep]
    
    process_grades, data, assign=assign, allassign=allassign, $
      weight=weight, class=class, semester=semester, test=test, $
      sendit=sendit, alldata=alldata, final=final, droplowest=droplowest
    srt = sort(alldata.current_grade)
    struct_print, alldata[srt]

return
end
