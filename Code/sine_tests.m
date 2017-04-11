A= -100+200.*rand(100000,1); % Initialize variables
arg0=pi/4;
arg1=sqrt(2);

t1=tic;
cos(A);
toc(t1)

t2=tic;
sin(A);
toc(t2)

t3=tic;
sin(A)+cos(A);
toc(t3)

t4=tic;
arg1*sin(A+arg0);
toc(t4)

t5=tic;
arg1*cos(A-arg0);
toc(t5)
