plot "barycenter_speed" u 2:3 w lines t "barycenter speed", \
"barycenter_coorinates" u 5:6:7 w lines t "moon", \
"barycenter_coorinates" u 8:9:10 w points t "sun", \
"barycenter_coorinates" u 11:12:13 w points t "barycenter"
pause -1
# "barycenter_coorinates" u ($8-$11):($9-$12):($10-$13) w lines t "sun-barycenter"
