name=["LF_CFL0.8","LF_CFL1.6","Leapfrog_CFL0.8"];
for k=1:3
    subplot(1,3,k);
    res=load(strcat(name(k),".csv"));
    plot(res(:,1),res(:,2),'o');
    hold on;
    % true res
    x=-2:0.01:3;
    len=length(x);
    y=zeros(len);
    for i=1:len
        if(abs(x(i)-0.8)<=1)
            y(i)=1-abs(x(i)-0.8);
        end
    end
    plot(x,y);
    title(strrep(name(k),'_',' ')); % String Replace _ for space
end
