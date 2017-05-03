clear all
x=importdata('prototype.txt');
y=[];
for i=1:length(x)
    if(x(i)~=0)
        y=[y;x(i)];
    end
end
histogram(y,200)
