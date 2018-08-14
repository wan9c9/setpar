prd=92
a_1=.661
b_1=-.388
a_2=.777
b_2=.333

a1=.5937
b1=-.2847
a2=.7042
b2=.2628


f(x)=a_1*cos(2.0*pi/prd*x)+b_1*sin(2.0*pi/prd*x)+a_2*cos(2.0*pi/prd*2.0*x)+b_2*sin(2.0*pi/prd*2.0*x)

g(x)=a1*cos(2.0*pi/prd*x)+b1*sin(2.0*pi/prd*x)+a2*cos(2.0*pi/prd*2.0*x)+b2*sin(2.0*pi/prd*2.0*x)
  
set title "Trend function"
set xrange [0:92]
#set yrange [-1:1]            
set size 1.0, 0.625
set terminal postscript portrait enhanced color dashed "Helvetica" 14 
set output "my-plot.ps"
plot f(x) title "Trend fitted by TPoiAR with trend", g(x) title "Trend fitted by PoiAR with trend"

#set terminal x11
#set size 1,1

