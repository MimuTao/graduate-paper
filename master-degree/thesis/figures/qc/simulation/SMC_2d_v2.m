%2D,
clc
clear
close all

%% System description.?
k1 = 2;    %mode number of system.
k2 = 2;    %mode number of controller.
final = 31;    %step
nh = 1;  %system dimension, nh->horitional nv->vertical.
nv = 1; 
n = nh+nv;
m = 2;  %input dimension
nu = 2;
ny = 2;


%% System Matrix
% feasible
% A{1} = [ -1.0 0.4; 0.2 -1.0];    A{2} = [ -0.8 0.6; 0.2 -1.2];   %A{2}=A{1}; 
% B{1} = [ 0.2 0.1; -0.2 -0.3];           B{2} = [ 0.2 -0.2; -0.1 -0.2];  %B{2}=B{1};  
% C{1} = [0.1 -0.1; -0.2 0.1];            C{2} = [-0.2 -0.1; 0.1 0.2];  
% D{1} = [-0.1 0; 0 0.1];              D{2} = [-0.1 0; 0 0.2];

A{1} = [ -1.0 0.6; -0.2 -1.3];    A{2} = [ -1 0.6; -0.2 -1.8];   %A{2}=A{1}; 
B{1} = [ 0.2 0.1; -0.2 -0.3];           B{2} = [ 0.2 -0.2; -0.1 -0.2];  %B{2}=B{1};  
C{1} = [1 0; 1.2 0.8];            C{2} = [1 0; 0.8 1];  
D{1} = [-0.1 0; 0 0.1];              D{2} = [-0.1 0; 0 0.2];


E{1} = [-0.1; 0.1];               E{2} = [-0.1; 0.1];

Z1 = 0.1;
Z2 = 0.1;

% A{1}=[0.6 0.4;0.4 0.5];               A{2}=[0.5 0.6; 0.5 0.3];  
% B{1}=[0.3 0.3;0.5 0.2];                          B{2}=[0.3 0.2;0.4 0.3];
% C{1}=[1 0;1 0.6];                     C{2}=[0.5 0.5; 0.8 0.6];        %C{2}=C{1};        
% D{1}=[0;0.3];                         D{2} = [0.1;0.2];     %  D{2}=D{1};   
% E{1}=[0.1;-0.1];                          E{2}=[0;-0.1];         

%% Mode Transfer Matrix.
P1 = [0.8 0.2; 0.3 0.7];    %system mode transfer function.
P2 = [0.6 0.4; 0.4 0.6];    %controller conditional mode transfer function.

%% LMI variables definition.
setlmis([]);
gamma_1 = lmivar(1,[1 0]); %represent the lambda max of Z^{T}R_{p}Z 
gamma_2 = lmivar(1,[1 0]);

% Q_ps T_ps
for p=1:k1
    for s=1:k2
        Q_rv{p,s} = lmivar(1,[n,1]);
    end
end
% R_p
T1 = [1;0];
T2 = [0;1];
for p=1:k1
    Rh_rv{p} = lmivar(1,[nh 1]);
    Rv_rv{p} = lmivar(1,[nv 1]);
%     R_rv{p} = lmivar(1,[nh 1;nv 1]);    %R_p = R^{-1}_p  
end
%Ks
for s=1:k2
    K_v{s} = lmivar(2,[nu n]);       
    L{s} = lmivar(2,[n n]);
end


%% Linear Matrix Inequality.

%condition1.
k=1;
for p=1:k1
    lmiterm([k 1 1 gamma_1],-1,1);
    lmiterm([k 1 2 0],Z1');
    lmiterm([k 2 2 Rh_rv{p}],-1,1);
    k = k+1;
end
for p=1:k1
    lmiterm([k 1 1 gamma_2],-1,1);
    lmiterm([k 1 2 0],Z2');
    lmiterm([k 2 2 Rv_rv{p}],-1,1);
    k = k+1;
end

%condition2.
for p=1:k1
    lmiterm([k 1 1 Rh_rv{p}],-T1,T1');
    lmiterm([k 1 1 Rv_rv{p}],-T2,T2');
    for s=1:k2
        lmiterm([k 1 1+s Rh_rv{p}],sqrt(P2(p,s))*T1,T1');
        lmiterm([k 1 1+s Rv_rv{p}],sqrt(P2(p,s))*T2,T2');
        lmiterm([k 1+s 1+s Q_rv{p,s}],-1,1);
    end
    k = k+1;
end

%condition3.
for p=1:k1
    for s=1:k2
        lmiterm([k 1 1 Q_rv{p,s}],1,1);
        lmiterm([k 1 1 L{s}],-1,1,'s');
        lmiterm([k 1 2 L{s}'],1,C{p}');
        lmiterm([k 1 2 K_v{s}'],1,D{p}');
        lmiterm([k 2 2 0],-1);
        for q=1:k1
            lmiterm([k 1 2+q L{s}'],sqrt(P1(p,q)),A{p}');
            lmiterm([k 1 2+q K_v{s}'],sqrt(P1(p,q)),B{p}');
            lmiterm([k 2+q 2+q Rh_rv{q}],-T1,T1');
            lmiterm([k 2+q 2+q Rv_rv{q}],-T2,T2');
        end
        k = k+1;
    end
end

lmisys = getlmis;

%% Parameter solution.
k_1 = 10;
k_2 = 10;
var_count = decnbr(lmisys);
c = [k_1;k_2;zeros(var_count-2,1)];
[xopt,bpot] = mincx(lmisys,c,[0 200 0 0 0 ]);
% [xopt,bpot] = mincx(lmisys,c,[10^(-6),1000,-1,200,0]);
result_gamma1 = dec2mat(lmisys,bpot,gamma_1);
result_gamma2 = dec2mat(lmisys,bpot,gamma_2);
gamma = k_1*result_gamma1 + k_2*result_gamma2
for s=1:k2
    Ks(:,:,s)= dec2mat(lmisys,bpot,K_v{s})*inv(dec2mat(lmisys,bpot,L{s}));
end


%% System mode
system_mode =[];
controller_mode = [];
system_mode(1,1) = 1;

%system mode transfer.
for j=1:final
    for i=1:final
        if(i*j~=1)
            if i==1
                flag1 = system_mode(i,j-1);
            else
                flag1 =system_mode(i-1,j);
            end
            if rand<P1(flag1,1)
                flag2 = 1;
            else 
                flag2 = 2;
            end
            system_mode(i,j) = flag2;
        end
    end
end

%controller mode transfer.
for i=1:final
    for j=1:final
        flag1 = system_mode(i,j);
        if rand<P2(flag1,1);
            flag2 = 1;
        else 
            flag2 = 2;
        end
        controller_mode(i,j) = flag2;
    end
end


% %% external disturbance
% w = [];
% for i=1:final
%     for j=1:final
%         if(i<=10)&&(j<=10)
%             w(i,j) = 0.2;
%         else
%             w(i,j) = 0;
%         end
%     end
% end
% 
%% system initial value
% x(h,//v,i,j)
x = [];
for j = 1:final
    if(j>=1)&&(j<=10)
        x(1,1,j) = 0.1;
    else
        x(1,1,j) = 0;
    end
end

for i = 1:final
    if(i>=1)&&(i<=10)
        x(2,i,1) = -0.1;
    else
        x(2,i,1) = 0;
    end
end

%% u = 0 时候的系统仿真
XX = [];
y = [];

for j = 1:final
    for i = 1:final
        flag3 = system_mode(i,j);
        AA = A{flag3};
        XX = AA*x(:,i,j); 
        x(1,i+1,j) = XX(1);
        x(2,i,j+1) = XX(2);
    end
end

% 3d-simulation plot
X1=[];             % horizontal state
X2=[];             % vertical state
for i=1:final
    for j=1:final
        X1(i,j)=x(1,i,j);
        X2(i,j)=x(2,i,j);
    end
end
x_m=linspace(0, final-1, final); % 在x轴上取final个点
y_m=linspace(0, final-1, final); % 在y轴上取final个点
[xx,yy]=meshgrid(x_m, y_m); % xx和yy都是finalxfinal的矩阵
mesh(xx,yy,X1);
xlabel('j')
ylabel('i')
zlabel('horizontal state without control input')

figure
mesh(xx,yy,X2);
xlabel('j')
ylabel('i')
zlabel('vertical state without control input')

%% u=Kx时，系统仿真
XX=[];
y=[];
for j=1:final
    for i=1:final
        flag3 = system_mode(i,j);
        AA = A{flag3};
        BB = B{flag3};
        CC = C{flag3};
        DD = D{flag3};
        EE = E{flag3};
        d1 = 0.3*norm(x(:,i,j));
        d = [d1;d1];
        
        
        flag4=controller_mode(i,j);
        KK=Ks(:,:,flag4);
        u(:,i,j)=KK*x(:,i,j);
        XX=AA*x(:,i,j)+BB*u(:,i,j); 
        x(1,i+1,j)=XX(1);
        x(2,i,j+1)=XX(2);
    end
end



% 画仿真三维图
X1=[];
X2=[];
U1=[];
U2=[];
for i=1:final
    for j=1:final
        X1(i,j)=x(1,i,j);
        X2(i,j)=x(2,i,j);
        U1(i,j)=u(1,i,j);
        U2(i,j)=u(2,i,j);
    end
end

x_m=linspace(0, final-1, final); % 在x轴上取101点
y_m=linspace(0, final-1, final); % 在y轴上取101点
[xx,yy]=meshgrid(x_m, y_m); % xx和yy都是101x101的矩阵

figure
mesh(xx,yy,X1);
xlabel('j')
ylabel('i')
zlabel('horizontal state with control input')

figure
mesh(xx,yy,X2);
xlabel('j')
ylabel('i')
zlabel('vertical state with control input')

figure
mesh(xx,yy,U1);
xlabel('j')
ylabel('i')
zlabel('horitical control input')
        
figure
mesh(xx,yy,U2);
xlabel('j')
ylabel('i')
zlabel('vertical control input')




% figure
% mesh(xx,yy,S1);
% xlabel('j')
% ylabel('i')
% zlabel('horitical sliding surface function s^{h}(i,j)')

% figure
% mesh(xx,yy,S2);
% xlabel('j')
% ylabel('i')
% zlabel('vertical sliding surface function s^{v}(i,j)')
% 
% figure
% grid
% for i=1:31
%     for j=1:31
%         if(system_mode(i,j)==1)
%             plot(i,j,'r+');hold on;
%         else
%             plot(i,j,'bx'); hold on;
%         end
%     end
% end
% xlabel('i')
% ylabel('j')
% legend('\itr(i, j) = 1','\itr(i, j) = 2');
% axis([0 32 0 32]);
% 
% figure
% for i=1:31
%     for j=1:31
%         if(controller_mode(i,j)==1)
%             plot(i,j,'bx');hold on;
%         else
%             plot(i,j,'r+'); hold on;
%         end
%     end
% end
% legend('\it\sigma(i, j) = 1','\it\sigma(i, j) = 2');
% % h = legend('\sigma(i,j) = 1','\sigma(i,j) = 2');
% % h = legend('sin(x)__','2*sin(x)__');
% % h1 = findobj(get(h,'Children'),'String','sin(x)__');
% % set(h1,'String','$sin(\hat{x})$','Interpreter','LaTex');
% xlabel('i')
% ylabel('j')
% axis([0 32 0 32]);
% grid on
%         
%         
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    

