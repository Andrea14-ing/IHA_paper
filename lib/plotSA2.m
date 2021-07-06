function h = plotSA2(SA,scale,col,varargin)
% Plots a set of Screw Axes into existing figure
% Input: SA = [nF x 6] pairs of n and s vectors of the axes
% plotSA2 : plot isas at steps of 10


nF = size(SA,1);
nF2 = size(scale,1); % if scale is provided as an array, pick the sample-wise value

if (nF2 ~= nF && nF2 ~= 1)
    error(' Inputs must be of the same length ');
end

if nF2 == 1
    for k=1:5:nF   
        n = SA(k,1:3);
        s = SA(k,4:6);
        P1 = s - scale.*n;
        P2 = s + scale.*n;
        h = plot3([P1(1);P2(1)],[P1(2);P2(2)],[P1(3);P2(3)] ,col,varargin{:});
        
    end
else
    for k=1:5:nF   
        n = SA(k,1:3);
        s = SA(k,4:6);
        P1 = s - scale(k).*n;
        P2 = s + scale(k).*n;
        h = plot3([P1(1);P2(1)],[P1(2);P2(2)],[P1(3);P2(3)] ,col,varargin{:});
        
    end
    
end

end
