function handle = EZquiver3(startVec, deltaVec, fig)

if nargin ==3
    figure(fig);
else
    figure()
end
hold on;
handle = quiver3(startVec(1),startVec(2),startVec(3), ...
                 deltaVec(1), deltaVec(2), deltaVec(3), ...
                 0);

end