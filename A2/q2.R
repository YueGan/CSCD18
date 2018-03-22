library(pracma)
gaze = c(3,-2,3)
w = gaze / sqrt(sum(gaze)^2)
w
t = c(0,0,1)
txc = cross(t,w)
u = txc / sqrt(sum(txc)^2)
u
v = cross(w,u)
v

Mcw = matrix(t(u), 3,1)
Mcw = cbind(Mcw, v)
Mcw = cbind(Mcw, w)
Mcw
pc1 = c(0.027836, -0.309684, 1.0)
pw = Mcw %*% (-7.0191 * pc1) + c(1,4,1)
pw

p2 = solve(Mcw) %*% (pw - c(1,6,1))
p2 = p2 / p2[3,1]
p2
