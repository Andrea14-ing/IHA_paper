function y = LogisticFunc(x,L,x0,k)
% Calculate the output of the logistic function for x [0,1]
y = L./(1 + exp(-k*(x-x0)));
end
