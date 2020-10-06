clear
global Pot_Mat len A_Vect;
num_trials = 10;

MEAobj = MCRackRecording('E:\Itamar\project\Biological Project\Recordings\Recording1.mcd');
T=MEAobj.getTrigger;  
M=MEAobj.getData(1:120,T{1}+200,30);
load layout_100_12x12.mat
% figure
% WaveformPhysicalSpacePlot(1:120,squeeze(M(:,1,:)),En,0)
mM=min(M,[],3);
% figure
% IntensityPhysicalSpacePlot(1:120,mM(:,1),En)
Pot_Mat = zeros(120);
for idx = 1:num_trials
    Pot_Mat = Pot_Mat + mM(:,120*(idx-1)+1:120*idx);
end
Pot_Mat = abs(Pot_Mat).*~eye(120)/num_trials;

%% remove bad electrodes - in this case 3, 64, 98, 119
bad_electrodes = [3 64 98 119];

for idx = sort(bad_electrodes,'descend')
    Pot_Mat(idx,:) = [];
    Pot_Mat(:,idx) = [];
end

%verify_matrix(Pot_Mat)

%% defining hyperparameters
dX=         200; %distance between electrodes
alpha0=     1.5; %exponential decay constant
noiseStd=   0.1; 
len =       length(Pot_Mat); %number of electrodes

%% Regression analysis of potentials and distances to see if data conforms to model
load('Solution_Coordinates.mat')
Coordinates = Coordinates *1e3;
for idx = sort(bad_electrodes,'descend')
    Coordinates(idx,:) =[];
end
r = dist(Coordinates');
v=Pot_Mat;
r(logical(eye(length(r)))) = 0;
v(logical(eye(length(Pot_Mat)))) = inf;


modelfun = @(c,x)(c(1)*(x.^-c(2)))+c(3);
modelfit = cell(len,1);
reg_Coef = zeros(len,3);
R_sq = zeros(len,1);
for idx = 1:len
    x = r(:,idx);
    x(idx) = [];
    y = v(:,idx);
    y(idx) = [];
    
    modelfit{idx} = fitnlm(x,y,modelfun,[0.5 0.5 0.5]);
    reg_Coef(idx,:) = modelfit{idx}.Coefficients.Estimate';
    R_sq(idx) = modelfit{idx}.Rsquared.Adjusted;
    
    %plot(x,modelfit{idx}.Fitted)
    
end
residuals = zeros(len-1);
for idx = 1:len
    tab = table2array(modelfit{idx}.Variables);
    [~,xorder]=sort(tab(:,1));
    restab = table2array(modelfit{idx}.Residuals(:,1));
    residuals(idx,:) = restab(xorder);
    
end
res_sum = [];
sum_size = [];
domain = [];
for idx = 1:len
   range = table2array(modelfit{idx}.Variables(:,1)); 
   res = table2array(modelfit{idx}.Residuals(:,1));
   for jdx = 1:len-1
       if isempty(find(domain==range(jdx)))
           dom_size = length(domain);
           domain(dom_size+1) = range(jdx);
           sum_size(dom_size+1) = 1;
           res_sum(dom_size+1) = res(jdx);
       else
           sum_size(find(domain==range(jdx))) = sum_size(find(domain==range(jdx))) +1;
           res_sum(find(domain==range(jdx))) = res_sum(find(domain==range(jdx))) + abs(res(jdx));
       end
   end
end
res_avg = res_sum./sum_size;
[~,domain_order] = sort(domain);
figure
plot(domain(domain_order),res_avg(domain_order))
           
    
plot(mean(abs(residuals)))
%%plotting a single regression
reg_num = 32;
figure
tab = table2array(modelfit{reg_num}.Variables);
[x,xorder]=sort(tab(:,1));
y = tab(xorder,2);
yhat = modelfit{reg_num}.Fitted(xorder);
plot(x,y,'.r')
hold on
plot(x,yhat)
title(['Regression of potential decay from electrode ' num2str(reg_num)])
xlabel('Distance [mm]')
ylabel('Potential [mV]')
%%
A =         [];
B =         [];
Aeq =       [];
beq =       [];
lb =        [-5000.*ones(len,3) , 1e-5*ones(len,1) ; 0 , -1000 , 0 , 0];
ub =        [5000.*ones(len,3) , 10000*ones(len,1) ; 3 , 1000 , 0 , 0];

nonlcon =   [];
options1 =  optimoptions('fmincon','Display','final-detailed','OptimalityTolerance',1e-20,'StepTolerance',1e-20,'ConstraintTolerance',1e-16,'MaxIterations',10000000,'MaxFunctionEvaluations',1000000,'Algorithm','sqp');
options2 =  optimoptions('fmincon','Display','final-detailed','OptimalityTolerance',1e-20,'StepTolerance',1e-20,'ConstraintTolerance',1e-16,'MaxIterations',100000000,'MaxFunctionEvaluations',100000000);
options3 =  optimoptions('fminunc','Display','final-detailed','OptimalityTolerance',1e-20,'MaxIterations',1000000,'MaxFunctionEvaluations',1000000);
options4 =  optimset('Display','final-detailed','MaxIter',1000000,'MaxFunEvals',100000,'TolFun',1e-10);


%%
% We will define a mesh grid on the z=0 plane as the initial conditions of
% the optimization. If the number of electrodes is not a perfect square we
% will create a mesh grid that with a larger number of electrodes and then
% % eliminate excess electrodes at random.
% Mesh_Size = nearest_upper_square(len);
% [X,Y]=meshgrid(0:dX:(dX*sqrt(Mesh_Size)-1),0:dX:(dX*sqrt(Mesh_Size)-1));
% Mesh=[X(:) Y(:) zeros(Mesh_Size,1) ones(Mesh_Size,1)]; %the mesh guess for the A constants is the mean of the constants calculated in the regression
% while length(Mesh) > len
%     Mesh(round(randinterval(1,1,1,1,120)),:) = [];
% end
% Mesh(:,4) = mean(reg_Coef(:,1));
% Mesh(end+1,:) =         [mean(reg_Coef(:,2)) 0 0 0]; 
% Pot_Mat = abs(Pot_Mat - mean(reg_Coef(:,3))).*~eye(size(Pot_Mat));
% A_struct = load('As_solved_by_optimization');
% A_vect = A_struct.Newt.Vector;
v0 = Pot_Mat;
v0(logical(eye(len))) = inf;
v0 = min(min(v0))-eps;
A_Vect = reg_Coef(:,1);
alpha = mean(reg_Coef(:,2));
near_sol = [Coordinates A_Vect ;alpha v0 0 0];
%% importing impedances from optimization and testing value at initial conditions with approximated impedances
Con_As = load('Con_As.mat');
GA_As = load('GA_As.mat');
Newt_As = load('Newt_As.mat');
PSO_As = load('PSO_As.mat');
SQP_As = load('SQP_As.mat');

Con_As.init_val = errorFun_3d([Coordinates Con_As.Con.Vector; alpha v0 0 0])
GA_As.init_val = errorFun_3d([Coordinates GA_As.GA.Vector'; alpha v0 0 0])
Newt_As.init_val = errorFun_3d([Coordinates Newt_As.Newt.Vector; alpha v0 0 0])
PSO_As.init_val = errorFun_3d([Coordinates PSO_As.PSO.Vector'; alpha v0 0 0])
SQP_As.init_val = errorFun_3d([Coordinates SQP_As.SQP.Vector; alpha v0 0 0])
regression_As.init_val = errorFun_3d([Coordinates A_Vect; alpha v0 0 0])
% x0 = [Coordinates A_Vect; alpha v0 0 0];
% Sol_Vector = fmincon(@(x) errorFun_3d(x),x0,A,B,Aeq,beq,lb,ub,nonlcon,options1);






%%
n=100;
nvars = 4*(len+1);
x0 =                zeros(len+1,4,n);
% Con.Vector =        zeros(len+1,4,n);
% Con.Time =          zeros(n,1);
% Con.Value =         zeros(n,1);
% Con.Dist =          zeros(n,1); 
% Con.Var_As =        zeros(n,1);
% Con.Mean_As =       zeros(n,1);

SQP.Vector =   zeros(len+1,4,n);
SQP.Time =     zeros(n,1);
SQP.Value =         zeros(n,1);
SQP.Dist =          zeros(n,1); 
SQP.Var_As =        zeros(n,1);
SQP.Mean_As =       zeros(n,1);

% Newt.Vector =   zeros(len+1,4,n);
% Newt.Time =     zeros(n,1);
% Newt.Value =         zeros(n,1);
% Newt.Dist =          zeros(n,1); 
% Newt.Var_As =        zeros(n,1);
% Newt.Mean_As =       zeros(n,1);
% 
% PSO.Vector =   zeros(len+1,4,n);
% PSO.Time =     zeros(n,1);
% PSO.Value =         zeros(n,1);
% PSO.Dist =          zeros(n,1); 
% PSO.Var_As =        zeros(n,1);
% PSO.Mean_As =       zeros(n,1);
% 
% GA.Vector =   zeros(len+1,4,n);
% GA.Time =     zeros(n,1);
% GA.Value =         zeros(n,1);
% GA.Dist =          zeros(n,1); 
% GA.Var_As =        zeros(n,1);
% GA.Mean_As =       zeros(n,1);


for idx = 1:n

    x0(:,:,idx) = near_sol+[noiseStd*randn(len,2) noiseStd*0.5*abs(randn(len,1)) noiseStd*abs(randn(len,1)) ; noiseStd*abs(randn) 0 0 0];
%     x0(:,:,idx) = near_sol;
    %Con(InteriorPoint)
%     tic
%     Con.Vector(:,:,idx) = fmincon(@(x) errorFun_3d(x),x0(:,:,idx),A,B,Aeq,beq,lb,ub,nonlcon,options2);
%     Con.Time(idx) = toc;
%     Con.Value(idx) = errorFun_3d(Con.Vector(:,:,idx));
%     Con.Dist(idx) = sum(sum(abs(dist(Con.Vector(1:len,1:3,idx)')-r)))/(-len+len^2);
%     Con.Var_As(idx) = var(Con.Vector(1:len,4,idx));
%     Con.Mean_As(idx) = mean(Con.Vector(1:len,4,idx));
%     
    %SQP
    tic
    SQP.Vector(:,:,idx) = fmincon(@(x) errorFun_3d(x),x0(:,:,idx),A,B,Aeq,beq,lb,ub,nonlcon,options1);
    SQP.Time(idx) = toc;
    SQP.Value(idx) = errorFun_3d(SQP.Vector(:,:,idx));
    SQP.Dist(idx) = sum(sum(abs(dist(SQP.Vector(1:len,1:3,idx)')-r)))/(-len+len^2);
    SQP.Var_As(idx) = var(SQP.Vector(1:len,4,idx));
    SQP.Mean_As(idx) = mean(SQP.Vector(1:len,4,idx));
    
    %Newt
%     tic
%     Newt.Vector(:,:,idx) = fminunc(@(x) errorFun_3d(x),x0(:,:,idx),options3);
%     Newt.Time(idx) = toc;
%     Newt.Value(idx) = errorFun_3d(Newt.Vector(:,:,idx));
%     Newt.Dist(idx) = sum(sum(abs(dist(Newt.Vector(1:len,1:3,idx)')-r)))/(-len+len^2);
%     Newt.Var_As(idx) = var(Newt.Vector(1:len,4,idx));
%     Newt.Mean_As(idx) = mean(Newt.Vector(1:len,4,idx));
     
%     %PSO
%     tic
%     PSO.Vector(:,:,idx) = vect2mat(particleswarm(@(x) errorFun_vector(x),nvars,mat2vect(lb),mat2vect(ub)));
%     PSO.Time(idx) = toc;
%     PSO.Value(idx) = errorFun_vector(mat2vect(PSO.Vector(:,:,idx)));
%     PSO.Dist(idx) = sum(sum(abs(dist(PSO.Vector(1:len,1:3,idx)')-r)))/(-len+len^2);
%     PSO.Var_As(idx) = var(PSO.Vector(1:len,4,idx));
%     PSO.Mean_As(idx) = mean(PSO.Vector(1:len,4,idx));
%      
%     %GA
%     tic
%     GA.Vector(:,:,idx) = vect2mat(ga(@(x) errorFun_vector(x),nvars,A,B,Aeq,beq,mat2vect(lb),mat2vect(ub)));
%     GA.Time(idx) = toc;
%     GA.Value(idx) = errorFun_vector(mat2vect(GA.Vector(:,:,idx)));
%     GA.Dist(idx) = sum(sum(abs(dist(GA.Vector(1:len,1:3,idx)')-r)))/(-len+len^2);
%     GA.Var_As(idx) = var(GA.Vector(1:len,4,idx));
%     GA.Mean_As(idx) = mean(GA.Vector(1:len,4,idx));
%     
     
    
    disp(['Now finishing iteration number:' num2str(idx)]) 
    
    
    
    
end
% [~,Con.best_iter] = min(Con.Value);
% [~,SQP.best_iter] = min(SQP.Value);
% [~,Newtbest_iter] = min(Newt.Value);
% [~,PSO.best_iter] = min(PSO.Value);
% [~,GA.best_iter] = min(GA.Value);

% 
% Best_Ans_Con = Con_Vector_iter(:,:,best_iter_Con);
% Best_Ans_SQP = SQP_Vector_iter(:,:,best_iter_SQP);
% Best_Ans_Newt = Newt_Vector_iter(:,:,best_iter_Newt);
% Best_Ans_PSO = PSO_Vector_iter(:,:,best_iter);
% Best_Ans_GA = GA_Vector_iter(:,:,best_iter_GA);

% 
% plot3(Best_Ans_Con(:,1),Best_Ans_Con(:,2),Best_Ans_Con(:,3),'.')
% plot3(Best_Ans_SQP(:,1),Best_Ans_SQP(:,2),Best_Ans_SQP(:,3),'.')
% plot3(Best_Ans_Newt(:,1),Best_Ans_Newt(:,2),Best_Ans_Newt(:,3),'.')
% plot3(Best_Ans_PSO(:,1),Best_Ans_PSO(:,2),Best_Ans_PSO(:,3),'.')
% plot3(Best_Ans_GA(:,1),Best_Ans_GA(:,2),Best_Ans_GA(:,3),'.')



% plot(1:length(Con_Value),sort(Con_Value))
% plot(1:n,Mean_As_Con)
% plot3(Con_Vector_iter(:,1,m),Con_Vector_iter(:,2,m),Con_Vector_iter(:,3,m),'.')

for idx = 1:10
    Near_Sol.f_x0(idx) = errorFun_3d(Near_Sol.x0(:,:,idx));
    
    
end