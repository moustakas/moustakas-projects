pro phys140_s14_grades, alldata, test=test, sendit=sendit, final=final
; jm14feb02siena - parse the grades for this class 

    path = getenv('TEACHING_DIR')+'/Phys140/140-S14/grades/'
    
    date = '14feb27' ; update this
    semester = 'Spring 2014'
    class = 'Physics 140 - General Physics II'

; read the grade spreadsheet downloaded from GoogleDocs
    gradefile = path+'phys140_grades_'+date+'.csv'
    data = read_gradefile(gradefile,unique_assignments=assign)

; specify the complete list of *possible* assignments and their
; relative weights
    allassign = ['Homework','Lab','Quizzes','Final Exam']
    weight = [0.25,0.15,0.40,0.20]
;   droplowest = [0,0,0,0]
    droplowest = [0,0,1,0]

;   keep = where(strmatch(data.last_name,'*Ippo*'))
;   keep = where(data.final_exam gt 0.0)
;   data = data[keep]
    
    process_grades, data, assign=assign, allassign=allassign, $
      weight=weight, class=class, semester=semester, test=test, $
      sendit=sendit, alldata=alldata, final=final, droplowest=droplowest
    srt = sort(alldata.current_grade)
    struct_print, alldata[srt]

return
end
