analit = open("analit.dat",'r')
method = open("trajectory",'r')
try:
    while(1):
        x = analit.readline()
        y = method.readline()
        x = list(map(float, x.split()))
        y = list(map(float, y.split()))
        print(x[0], end=" ")
        diff = -1 * (x[1] ** 2 + x[2] ** 2 + x[3] ** 2 - y[1] ** 2 - y[2] ** 2 - y[3] ** 2) / (y[1] ** 2 + y[2] ** 2 + y[3] ** 2)
        print(diff)
finally:
    analit.close()
    method.close()
