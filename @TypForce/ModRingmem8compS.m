function [S] = ModRingmem8compS(r1,r2,z1,z2)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('r1', @(x) isnumeric(x));
ip.addRequired('r2', @(x) isnumeric(x));
ip.addRequired('z1', @(x) isnumeric(x));
ip.addRequired('z2', @(x) isnumeric(x));
ip.parse(r1,r2,z1,z2);
%----------------------------------------------------------------------------------------
    if r1>r2
        r_tem=r1;
        r1=r2;
        r2=r_tem;
    end
    dz=abs(z1-z2);
    S=pi*(r2*sqrt(r2^2+dz^2/(1-r1/r2)^2)-r1*sqrt(r1^2+(r1^2/r2^2)*dz^2/(1-r1/r2)^2));
end


