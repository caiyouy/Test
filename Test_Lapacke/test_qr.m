k=1e7;
A=[k/sqrt(2)+0.5,-k/sqrt(6)+sqrt(3)/2;...
   -k/sqrt(2),k/sqrt(6);...
   k/sqrt(2)-0.5,-k/sqrt(6)-sqrt(3)/2];
AA=A'*A

b=[3;2;-1];
Ab=A'*b