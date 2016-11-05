pro phys310_f16_grades, alldata, test=test, sendit=sendit, final=final
; jm15sep20 - parse the grades for this class 

    path = getenv('TEACHING_DIR')+'/310-F16/grades/'
    
    date = '16oct30' ; update this
    semester = 'Fall 2016'
    class = 'PHYS310 - Mechanics I'

; read the grade spreadsheet downloaded from GoogleDocs
    gradefile = path+'phys310_grades_'+date+'.csv'
    data = read_gradefile(gradefile,unique_assignments=assign)

; specify the complete list of *possible* assignments and their
; relative weights
    allassign = ['Homework','Laboratory','Quizzes','Final']
    weight = [0.25,0.25,0.30,0.20]
;   droplowest = [1,0,1,0]
    droplowest = [1,0,1,0]

;   data.quizzes *= 1.05
;   data.final *= 1.03

;   keep = where(data.final_exam gt 0.0)
    keep = where(strmatch(data.first_name,'*Patrick*'))
    data = data[keep]

    process_grades, data, assign=assign, allassign=allassign, $
      weight=weight, class=class, semester=semester, test=test, $
      sendit=sendit, alldata=alldata, final=final, droplowest=droplowest
    srt = sort(alldata.current_grade)
    struct_print, alldata[srt]

return
end
