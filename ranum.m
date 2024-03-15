function rnums=ranum(rprime,mean,sigm,m,n)
% Use following prime numbers
% .8191,.2047,.131071,.54287
% Enter the mean and standard deviation of the every row of the matrix.

rnums=zeros(m,n);
for j=1:m   
    for i=1:n
        rq=urand(rprime);
        rs=urand(rq);
        rnum=sigm(m,m)*sqrt(-2*log(rq))*cos(6.28*rs)+mean;
        rnums(j,i)=rnum;
        rprime=rs;
    end
end
end



%%
function unum=urand(g)
d=331*g;
i=floor(d);
unum=d-i;
end
%%