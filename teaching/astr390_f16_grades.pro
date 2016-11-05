pro astr390_f16_grades, alldata, test=test, sendit=sendit, final=final
; jm15jan31siena - parse the grades for this class 

    path = getenv('TEACHING_DIR')+'/390-F16/grades/'
    
    date = '16oct24' ; update this
    semester = 'Fall 2016'
    class = 'ASTR390 - Principles of Astrophysics I'

; read the grade spreadsheet downloaded from GoogleDocs
    gradefile = path+'astr390_grades_'+date+'.csv'
    data = read_gradefile(gradefile,unique_assignments=assign)

; specify the complete list of *possible* assignments and their
; relative weights
    allassign = ['Homework','Quizzes','Talk','Final Problem Set']
    weight = [0.30,0.40,0.15,0.15]
    droplowest = [1,1,0,0]

    keep = where(strmatch(data.first_name,'*Alexander*'))
;   keep = where(strmatch(data.first_name,'*Cristian*'))
    data = data[keep]

    process_grades, data, assign=assign, allassign=allassign, $
      weight=weight, class=class, semester=semester, test=test, $
      sendit=sendit, alldata=alldata, final=final, droplowest=droplowest
    srt = sort(alldata.current_grade)
    struct_print, alldata[srt]

return
end
