function [varargout] = nuclearSize(t,method,n,phase)
r=0; a=0; b=0;

L = 270; %half-embryo linear length in microns
h = 25; %compartment height, approx. in microns
switch phase
    case {'interphase','mitosis'}
        switch n
            case 13
                NC = 10;
                Tspani = 0; Tspanf = 3.6072;
            case 19
                NC = 11;
                Tspani = 7.69538; Tspanf = 11.5431;
            case 26
                NC = 12;
                Tspani = 15.6312; Tspanf = 21.8838;
            case 36
                NC = 13;
                Tspani = 25.9719; Tspanf = 37.9960;
            case 51
                NC = 14;
                Tspani = 43.2953; Tspanf =  97.6954;
        end
    case {'interphase2D','mitosis2D'}
        switch n{1}
            case {5,'NC10','NC10m'}
                n=n{2};
                NC = 10;
                Tspani = 0; Tspanf = 3.6072;
            case {8,'NC11','NC11m'}
                n=n{2};
                NC = 11;
                Tspani = 7.69538; Tspanf = 11.5431;
            case {12,'NC12','NC12m'}
                n=n{2};
                NC = 12;
                Tspani = 15.6312; Tspanf = 21.8838;
            case {17,'NC13','NC13m'}
                n=n{2};
                NC = 13;
                Tspani = 25.9719; Tspanf = 37.9960;
            case {25,'NC14'}
                n=n{2};
                NC = 14;
                Tspani = 43.2953; Tspanf =  97.6954;
        end

end
w = L/n; %width of a compartment
Vol = w*w*h; %volume of a compartment
Am = 2*w*h; % compartment area available for transport
Acs = w*w; % Area available for Toll receptors

% scaled time
if strcmp(phase,'premitosis')
    t_hat =1;
    NC = 10;
    Tspani = 0; Tspanf = 3.6072;
else
    t_hat = (t-Tspani)/(Tspanf-Tspani);
end
% area and volume of a nucleus in squred or cubed microns
switch method
    case {'dynamic','dynamic2D'}
        switch NC
            case 10
                %NC10
                r = 1.8246*t_hat + 3.3484;
                An = 4*pi*r^2;
                Vn = 4/3*pi*r^3;
                
                drdt = 1.8;
                dVndt = 4*pi*r^2*drdt;
            case 11
                %NC11
                r = 2.2070*t_hat + 3.2002;
                An = 4*pi*r^2;
                Vn = 4/3*pi*r^3;
                
                drdt = 2.2;
                dVndt = 4*pi*r^2*drdt;
                
            case 12
                %NC12
                r = 1.7224*t_hat + 2.7659;
                An = 4*pi*r^2;
                Vn = 4/3*pi*r^3;
                
                drdt = -2.4*t_hat+3;
                dVndt = 4*pi*r^2*drdt;
            case 13
                % NC13
                y0 = 1.9070;
                V = 2.9922;
                K = 0.1934;
                r = y0+V*t_hat./(K+t_hat);
                drdt = (K*V)/(K+t_hat)^2;
%                 
%                 4th order polynomial does cause high error
%                 r = -9.3*t_hat^4+24*t_hat^3-23*t_hat^2+11*t_hat+2;
%                 drdt = -37.2*t_hat^3+72*t_hat^2-46*t_hat+11;
%                 
%                 3rd order polynomial
% r = 5.4*t_hat^3-11*t_hat^2+8.3*t_hat+2.1;
% drdt = 16.2*t_hat^2-22*t_hat+8.3;
%
%                 
%                 quadratic does not cause high error
%                 r = -3.2*t_hat^2+5*t_hat+2.4;
%                 drdt = -6.4*t_hat+5;
%                 
%                 linear, does not cause high error
%                 r = 1.8*t_hat+2.9;
%                 drdt = 1.8;
                
                An = 4*pi*r^2;
                Vn = 4/3*pi*r^3;
                
                dVndt = 4*pi*r^2*drdt;
                
            case 14
                %NC14 - a spheroid with major and minor axes
                
                y0 = 1.6600;
                V = 1.3310;
                K = 0.0132;
                a = y0+V*t_hat./(K+t_hat);
                b = y0+V*t_hat./(K+t_hat)+2*t_hat; 
                % this is approximated from live imaging
%                 a = 0.29*t_hat+2.8;
%                 b = 2.3*t_hat+2.8;
                if a<b
                    e = sqrt(1-(a^2/b^2));
                    An = 2*pi*a^2*(1+b/a/e*asin(e));
                else
                    An = 4*pi*b^2;
                end
                Vn = 4/3*pi*a^2*b;
                dadt = K*V/(K+t_hat)^2;
                dbdt = K*V/(K+t_hat)^2+2;
% dadt = 0.29;
% dbdt = 2.3;
                dVndt = 4/3*pi*a*(2*b*dadt+a*dbdt);
                
        end
    case {'static','static2D'}
        switch NC
            case 10
                %NC10
                r = 4.2;
                An = 4*pi*r^2;
                Vn = 4/3*pi*r^3;
                dVndt = 0;
            case 11
                %NC11
                r = 4.3;
                An = 4*pi*r^2;
                Vn = 4/3*pi*r^3;
                dVndt = 0;
                
            case 12
                %NC12
                r = 3.6980;
                An = 4*pi*r^2;
                Vn = 4/3*pi*r^3;
                dVndt = 0;
            case 13
                % NC13
                r = 3.8391;
                An = 4*pi*r^2;
                Vn = 4/3*pi*r^3;
                dVndt = 0;
                
            case 14
                %NC14 - a spheroid with major and minor axes
                b = 2.9080;
                a = 3.9080; % this is approximated from live imaging
                
                if b<a
                    e = sqrt(1-(b^2/a^2));
                    p1 = atanh(e);
                    p2 = (1-e^2)/e;
                    p3 = 2*pi*a^2;
                    An = p3*(1+p2*p1);
                    
                else
                    An = 4*pi*b^2;
                end
                Vn = 4/3*pi*a^2*b;
                dVndt = 0;
                
        end
end

% volume of compartment - volume of nucleus
Vc = Vol - Vn;
if Vc < 0
    error('Negative Volume!')
end


% scaling parameters
% An14 = 231.8830;
% Vn14 = 308.1286;
An14 = 160.1618;
Vn14 = 186.034;
% scaled parameters
Vol = Vol/Vn14;
An = An/An14;
Am = Am/An14;
Acs = Acs/An14;
Vn = Vn/Vn14;
dVndt = dVndt/(Tspanf-Tspani)/Vn14;
% Vn = max(0.5,Vn/Vn14);
% dVndt = min(dVndt/Vn14,3);
Vc = Vc/Vn14;

if strcmp(phase,'mitosis')
    r=0;
    a=0;
    b=0;
    An = 0;
    dVndt = 0;
    Vn = 0;
    Vc = Vol;
end
% if Vn < 1e-2
%     lastwarn('Vn -> 0 inside interphase')
% end
argout = {An,Am,Acs,Vn,Vc,dVndt,r,a,b,Vol};
varargout = cell(1,nargout);

for i = 1:nargout
    varargout{i} = argout{i};
end

end

