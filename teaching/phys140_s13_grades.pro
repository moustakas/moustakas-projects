pro phys140_s13_grades, alldata, test=test, sendit=sendit, final=final
; jm12oct05siena - parse the grades for this class 

    path = getenv('TEACHING_DIR')+'/140-S13/grades/'
    
    date = '13may13' ; update this
    semester = 'Spring 2013'
    class = 'Physics 140 - General Physics II'

; read the grade spreadsheet downloaded from GoogleDocs
    gradefile = path+'phys140_grades_'+date+'.csv'
;   gradefile = path+'junk.csv'
    data = read_gradefile(gradefile,unique_assignments=assign)

; specify the complete list of *possible* assignments and their
; relative weights
    allassign = ['Homework','Lab','Midterm 1','Midterm 2','Midterm 3','Final Exam']
    weight = [0.2,0.15,0.15,0.15,0.15,0.2]

    keep = where(data.final_exam gt 0.0)
    data = data[keep]
    
    process_grades, data, assign=assign, allassign=allassign, $
      weight=weight, class=class, semester=semester, test=test, $
      sendit=sendit, alldata=alldata, final=final
    srt = sort(alldata.current_grade)
    struct_print, alldata[srt]

return
end
