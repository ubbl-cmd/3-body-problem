splot "barycenter_coorinates" u 2:3:4 w lines t "barycenter speed", \
"barycenter_coorinates" u 5:6:7 w lines t "moon", \
"barycenter_coorinates" u 8:9:10 w points t "sun", \
"barycenter_coorinates" u 11:12:13 w points t "barycenter"
pause -1
