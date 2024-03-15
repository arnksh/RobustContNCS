function [K,P]=dcontr122(a,bm,g)
n=length(a);
e1=0.0016;
e2=0.00018;
setlmis([])
p=lmivar(1,[n 1]);
phi=lmivar(2,[1 n]);

for i=1:2
    lmiterm([i 1 1 p],-e1,1);
    lmiterm([i 2 1 p],a*bm(i),1);
    lmiterm([i 2 1 phi],g,1);
    lmiterm([i 2 2 p],-1/e1,1);
    lmiterm([i 3 3 p],-e2,1);
    lmiterm([i 4 1 p],a,1);
    lmiterm([i 4 3 p],a*bm(i),1);
    lmiterm([i 4 3 phi],g,1);
    lmiterm([i 4 4 p],-1/e2,1);
end
    lmiterm([3 1 1 p],-1,1);
    
    lmis=getlmis;
    [~,x]=feasp(lmis);
    P=dec2mat(lmis,x,p);
    Phi=dec2mat(lmis,x,phi);
    K=Phi/P;
end
