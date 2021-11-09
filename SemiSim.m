% Simulating a wheel rolling on the ground
% Calculate moment of inertia tensor
M = 5;
N = 100;
r= 0.5;
g = 9.8;
S = 9999;
mu = 0.5;
Tf = 10;
thet = (0:N-1)'*2*pi/N;
Omega = -[0 0 1];
pos = [r*sin(thet) r*cos(thet) zeros(N,1)]+[0.5 1+r 0];
x_cm = mean(pos);
x_tilde = pos-x_cm;
u_cm = [0.5 0 0];
I = zeros(3);
for k=1:N
    I = I+M/N*(norm(x_tilde(k,:))^2*eye(3)-x_tilde(k,:)'*x_tilde(k,:));
end
% I does not change for a circle
L = I*Omega';
dt = 1e-2;
nSteps= Tf/dt;
xcms = zeros(nSteps+1,3);
xcms(1,:)=x_cm;
Omegas = zeros(nSteps,1);

for iStep=1:nSteps
% Compute I and angular velocity from L
% I = zeros(3);
% for k=1:N
%     I = I+M/N*(norm(x_tilde(k,:))^2*eye(3)-x_tilde(k,:)'*x_tilde(k,:));
% end 
% I
% Compute forces on every pt (need real xs for this)
% Numerical method
Omega = (I \ L)';
Omegas(iStep)= Omega(3);
pos = x_tilde+x_cm;
u = cross(Omega.*ones(N,3),x_tilde)+u_cm;
F = zeros(N,3);
for iPt=1:N
    F(iPt,:) = GroundForce3D(pos(iPt,:),u(iPt,:),mu,S,M*g/N);
end
TotalF = sum(F);

% Rotate the Xks
for iPt=1:N
    x_tilde(iPt,:)=rotate(x_tilde(iPt,:),Omega*dt);
end
% Update center of mass
u_cm = dt/M*TotalF+u_cm;
L = L+dt*sum(cross(x_tilde,F))';
x_cm = x_cm+dt*u_cm;
xcms(iStep+1,:)=x_cm;
plot([pos(:,1);pos(1,1)],[pos(:,2);pos(1,2)])
hold on
plot(pos(3,1),pos(3,2),'ro')

x1=1:0.01:2;
y1=-sqrt(1-(x1-2).*(x1-2))+1;
x2=2:0.1:3;
y2=0;
x3=3:0.01:4;
y3=-sqrt(1-(x3-3).*(x3-3))+1;
x4=0:0.1:1;
y4=1;
plot(x1,y1)
plot(x2,y2)
plot(x3,y3)
plot(x4,y4)

%plot(xlim,[0 0],'-k')
axis ([-1 7 0 8])
drawnow
hold off
end

function force = GroundForce3D(x,u,mu,S,mg)
    [hOfRamp,gradH] = Hvals3dSEMI(x);
    force = -mg*[0 1 0];
    if (hOfRamp <= 0)
        % Add the ground forces
        n = gradH/norm(gradH); 
        disp = -(hOfRamp)/norm(gradH);
        Utan = u - dot(u,n)*n;
        if (norm(Utan) > 1e-10)
            UtanHat = Utan/norm(Utan);
        else
            UtanHat = [0 0 0];
        end
        force = force + S*(disp*(n-mu*UtanHat));
    end
end
    
function [h,gradRamp] = Hvals3dSEMI(pos) % specifies the ramp
    func1 = -sqrt(1-(pos(1)-2)^2)+1;
    func2 = 0;
    func3 = -sqrt(1-(pos(1)-3)^2)+1;
    func4 = 1;
if(pos(1)>=3.7||2 <= pos(1) && pos(1)<= 3)
    h = pos(2)-func2;
    gradRamp = [0 1 0];
elseif(pos(1)<=1)
    h = pos(2)-func4;
    gradRamp = [0 1 0];
elseif(pos(1)<2&& 1 < pos(1))
    deriv = (pos(1)-2)*(1-(pos(1)-2)^2)^(-1/2);
    h = pos(2)-func1;
    gradRamp = [-deriv 1 0];
elseif(pos(1)>3&&pos(1)<3.7)
    deriv = (pos(1)-3)*(1-(pos(1)-3)^2)^(-1/2);;
    h = pos(2)-func3;
    gradRamp = [-deriv 1 0];
end  
end


% function [h,gradRamp,ramPlot] = Hvals3dSEMI(pos) % specifies the ramp
%     func1 = (pos(1)-1)^2;    
%     func2 = 0;
%     deriv = 2*(pos(1)-1);
%     h = func1;
%     gradRamp = [1 -1/deriv 0]/norm([1 -1/deriv 0]);
%     
%     syms x
%     ramPlot = piecewise(x>0, func1, func2);
% %     
% end


    
function rotated_x = rotate(x, Om)
    nOm = norm(Om);
    Omhat = Om/nOm;
    Px = Omhat*dot(Omhat,x);
    rotated_x = Px+cos(nOm)*(x-Px)+sin(nOm)*cross(Omhat,x);
    if(nOm < 1e-10)
        rotated_x = x;
    end
end
